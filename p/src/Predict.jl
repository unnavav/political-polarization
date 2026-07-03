module Predict

using Random
using Printf: @printf, @sprintf 
using Statistics: mean

using ..Solvers: KSsolver
using ..DistrTools: getDistr, transitDistr
using ..Compute: weight, supnorm

export simz, transition, perfectForesight, genForecastData, update_forecast, run_KS

function simz(N_T, Nz, rnseed, P_mat)
    Random.seed!(rnseed)
    efsim = rand(N_T)

    izsim = zeros(Int, N_T)
    izsim[1] = div(Nz, 2) + 1        # middle-ish index

	# here we normalize the transition matrix to ensure that it sums to 1 across rows
	P_mat = P_mat ./ sum(P_mat, dims=2) 

    cumP_mat = cumsum(P_mat, dims=2)

    for t = 1:N_T-1
        u = efsim[t+1]
        csumvec = @view cumP_mat[izsim[t], :]
        izsim[t+1] = findfirst(>=(u), csumvec)
    end

    return izsim
end

function genForecastData(V, V0, G, G0, C, Kfore, params, policies, prices, NT, rnseed, vTol)

	nk, nz, nl, na = size(V);
	Kgrid = params.Kgrid; 

	futureKs = exp.(Kfore[:, 1] .+ Kfore[:, 2] .* log.(Kgrid)')   # nz × nk
    @assert size(futureKs) == (nz, nk)

	CI_hh = CartesianIndices((1, nz))     # for solve to hack the current setup
    V, G, C, ~, ~ = KSsolver(V, V0, G, G0, C, futureKs, Kgrid, params, policies,
                       prices, CI_hh, vTol)

	zt = simz(NT, nz, rnseed, params.π_z)

	# choose a random starting point for the simulation
	ik0 = cld(nk, 2); K0 = Kgrid[ik0];
    Kt = zeros(NT+1); Kt[1] = Kgrid[ik0]

	CI_dist = CartesianIndices(params.π_z)   # (nz,nz) for getDistr
	LI_dist = LinearIndices(params.π_z)
    # initial distribution: stationary at starting K, lifted to pair-state
    G_start = G[ik0, :, :, :]                        # (nz,nl,na) single-z
	G_pair = zeros(nz*nz, nl, na)
	for it in 1:(nz*nz)
		today = CI_dist[it][2]
		G_pair[it, :, :] = G_start[today, :, :]
	end

    # ... lift to pair-state (nz²) if getDistr is pair-state ...
    μ_prev, _ = getDistr(G_pair, params.amu, params.agrid, params.π_l, params.π_z,
                         CI_dist, LI_dist, params.ϕ)

	#@printf("getDistr: sum(μ_prev) = %.10f  (want 1.0)\n", sum(μ_prev))
	#@printf("          min = %.3e  (want ≥ 0, no negatives)\n", minimum(μ_prev))

	# μ_pair is (nz², nl, nmu); collapse to (nz, nl, nmu) by today's z (all for when z transition matters)
	nmu = length(params.amu);
	μ_today = zeros(nz, nl, nmu)
	for it in 1:(nz*nz)
		today = CI_dist[it][2]
		μ_today[today, :, :] .+= μ_prev[it, :, :]
	end

	μ_transit = μ_today[zt[1]:zt[1], :, :];
	#@printf("collapse: sum(μ_today) = %.10f  (want = sum(μ_prev))\n", sum(μ_today))

	for t in 1:NT
		if t%2500 == 0
			@printf("\tSimulating period %i of %i\n", t, NT)
		end

		K = Kt[t]; ix, we = weight(Kgrid, K);
		iz = zt[t];  # today's TFP state
		G_t = we .* G[ix, iz, :, :] .+ (1-we) .* G[ix+1, iz, :, :];   # (nl, na)

		# update the distribution for the next period
		mass_before = sum(μ_transit)
		μ_transit, Kt[t+1] = transitDistr(G_t, μ_transit, params.amu, params.agrid, params.ϕ, params.π_l)
		mass_after = sum(μ_transit)

		if abs(mass_after - mass_before) > 1e-8
			@printf("LEAK at t=%i: before=%.10f after=%.10f  Δ=%.3e  (K=%.4f)\n",
					t, mass_before, mass_after, mass_after - mass_before, K)
		end

	end
    
	return Kt, zt

end

function update_forecast(Kt, zt, nz, burn_in)
    Kfore_new = zeros(nz, 2)
    R2 = fill(NaN, nz)
    counts = zeros(Int, nz)

    for z in 1:nz
        # periods where TODAY's state is z, post burn-in, with a valid t+1
        idx = [t for t in (burn_in+1):(length(Kt)-1) if zt[t] == z]
        counts[z] = length(idx)

        if length(idx) < 5          # too few to estimate a 2-param rule
            Kfore_new[z, :] = [0.0, 1.0]   # fallback: identity in logs
            continue
        end

        x = log.(Kt[idx])            # log K_t
        y = log.(Kt[idx .+ 1])       # log K_{t+1}
        X = hcat(ones(length(x)), x) # design matrix [1  logK]
        β = X \ y                    # OLS: [intercept, slope]
        Kfore_new[z, :] = β

        ŷ = X * β
        ss_res = sum((y .- ŷ).^2)
        ss_tot = sum((y .- mean(y)).^2)
        R2[z] = ss_tot > 0 ? 1 - ss_res/ss_tot : NaN
    end

    return Kfore_new, R2, counts
end

function run_KS(V, V0, G, G0, C, params, policies, prices,
                NT, rnseed, vTol, dTol; burnin=500, λ_damp=0.3, maxout=100)

    nk, nz, nl, na = size(V)
    Kfore = repeat([0.0 1.0], nz, 1)   # nz×2, start at log-identity
    foredist = 1e5
    outer_ct = 1

    while foredist > dTol && outer_ct ≤ maxout
        # solve HH + simulate under current Kfore
        Kt, zt = genForecastData(V, V0, G, G0, C, Kfore, params, policies,
                                 prices, NT, rnseed, vTol)

        Kfore_new, R2, counts = update_forecast(Kt, zt, nz, burnin)

		println("\nForecast rules:  log K' = a + b·log K")
		println("─"^58)
		@printf("  %-8s %11s %11s %9s %8s\n", "z-state", "a", "b", "R²", "n")
		println("─"^58)
		for z in 1:nz
			flag = R2[z] < 0.99 ? "  ⚠" : ""
			@printf("  %-8d %11.6f %11.6f %9.4f %8d%s\n",
				z, Kfore_new[z,1], Kfore_new[z,2], R2[z], counts[z], flag)
		end
		println("─"^58)

        foredist = maximum(abs.(Kfore_new .- Kfore))

        @printf("Outer %2i | foredist = %.6f | R² = [%s] | counts = [%s]\n",
                outer_ct, foredist,
                join([@sprintf("%.4f", r) for r in R2], ", "),
                join(string.(counts), ", "))

        # damped update
        Kfore = λ_damp .* Kfore_new .+ (1 - λ_damp) .* Kfore
        outer_ct += 1
    end

    if foredist ≤ dTol
        @printf("\nConverged in %i outer iters. foredist = %.6f\n", outer_ct-1, foredist)
    else
        @printf("\nHit maxout=%i without converging. foredist = %.6f\n", maxout, foredist)
    end

    return Kfore, Kt, zt
end


end # module