module Predict

using Random
using Printf: @printf, @sprintf 
using Statistics: mean, median

using ..Solvers: KSsolver
using ..DistrTools: getDistr, transitDistr
using ..Compute: weight, supnorm, summarizeKtByTransition

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



function genForecastData(V, V0, G, G0, C, Kfore, params, policies, prices, zt, vTol;
    verbose = false)

	CI = CartesianIndices(params.π_z)   # (nz,nz) for getDistr
	LI = LinearIndices(params.π_z)

	nk, nt, nl, na = size(V);
	Kgrid = params.Kgrid; 

	futureKs = exp.(Kfore[:, 1] .+ Kfore[:, 2] .* log.(Kgrid)')   # nt × nk
    @assert size(futureKs) == (nt, nk)

    V, G, C, ~, ~ = KSsolver(V, V0, G, G0, C, futureKs, Kgrid, params, policies,
                       prices, CI, LI, vTol; verbose)

	# choose a random starting point for the simulation
    NT = length(zt);
    ik0 = cld(nk, 2); K0 = Kgrid[ik0];
    Kt = zeros(NT+1); Kt[1] = Kgrid[ik0]

    # initial distribution: stationary at starting K, lifted to pair-state
    G_start = G[ik0, :, :, :]                        # (nt, nl,na) (base K)


    # ... lift to pair-state (nz²) if getDistr is pair-state ...
    μ_transit, _ = getDistr(G_start, params.amu, params.agrid, params.π_l, params.π_z,
                         CI, LI, params.ϕ)

	#@printf("getDistr: sum(μ_prev) = %.10f  (want 1.0)\n", sum(μ_prev))
	#@printf("          min = %.3e  (want ≥ 0, no negatives)\n", minimum(μ_prev))

	#@printf("collapse: sum(μ_today) = %.10f  (want = sum(μ_prev))\n", sum(μ_today))
    med_ind = ceil(Int, median(1:length(params.zgrid)))
    it_t = zeros(Int, NT);
    it_t[1] = LI[med_ind, med_ind]   # initial pair-state index
    it_t[2:NT] = [LI[zt[t-1], zt[t]] for t in 2:NT];
    μ_transit = μ_transit[it_t[1], :, :]

	for t in 1:NT
		if t%2500 == 0 && verbose
			@printf("\tSimulating period %i of %i\n", t, NT)
		end

		K = Kt[t]; ix, we = weight(Kgrid, K);
		it = it_t[t] # getting which z transition we're in
		G_t = we .* G[ix, it, :, :] .+ (1-we) .* G[ix+1, it, :, :];   # (nl, na)

		# update the distribution for the next period
		mass_before = sum(μ_transit)
		μ_transit, Kt[t+1] = transitDistr(G_t, μ_transit, params.amu, params.agrid, params.ϕ, params.π_l)
		mass_after = sum(μ_transit)

		if abs(mass_after - mass_before) > 1e-8
			@printf("LEAK at t=%i: before=%.10f after=%.10f  Δ=%.3e  (K=%.4f)\n",
					t, mass_before, mass_after, mass_after - mass_before, K)
		end

	end
    
	return Kt, it_t

end

function update_forecast(Kt, it_t, nt, burn_in)
    Kfore_new = zeros(Float64, nt, 2)
    R2 = fill(NaN, nt)
    counts = zeros(Int, nt)

    for z1z2 in 1:nt
        # periods where TODAY's state is z, post burn-in, with a valid t+1
        idx = [t for t in (burn_in+1):(length(Kt)-1) if it_t[t] == z1z2]
        counts[z1z2] = length(idx)

        if length(idx) < 5          # too few to estimate a 2-param rule
            Kfore_new[z1z2, :] = [0.0, 1.0]   # fallback: identity in logs
            continue
        end

        x = log.(Kt[idx])            # log K_t
        y = log.(Kt[idx .+ 1])       # log K_{t+1}
        X = hcat(ones(length(x)), x) # design matrix [1  logK]
        β = X \ y                    # OLS: [intercept, slope]
        Kfore_new[z1z2, :] = β

        ŷ = X * β
        ss_res = sum((y .- ŷ).^2)
        ss_tot = sum((y .- mean(y)).^2)
        R2[z1z2] = ss_tot > 0 ? 1 - ss_res/ss_tot : NaN
    end

    return Kfore_new, R2, counts
end

function run_KS(V, V0, G, G0, C, params, policies, prices,
                zt, Kfore, vTol, dTol; burnin=500, λ_damp=0.3, maxout=100,
                verbose = false)

    _, nt, _, _ = size(V)
    foredist = 10.0
    outer_ct = 1

    NT = length(zt);
    Kt = zeros(NT+1); 

    while foredist > dTol && outer_ct ≤ maxout

        # vTol_outer = max(vTol, foredist * 1e-2)
        vTol_outer = vTol;
        
        # solve HH + simulate under current Kfore
        Kt, it_t = genForecastData(V, V0, G, G0, C, Kfore, params, policies,
                                 prices, zt, vTol_outer, verbose = verbose)

        Kfore_new, R2, counts = update_forecast(Kt, it_t, nt, burnin)
        foredist = maximum(abs.(Kfore_new .- Kfore))

        if verbose
            CI = CartesianIndices(params.π_z)
            println("\nForecast rules:  log K' = a + b·log K")
            println("─"^62)
            @printf("  %-10s %11s %11s %9s %8s\n", "z₋₁→z", "a", "b", "R²", "n")
            println("─"^62)
            for z1z2 in 1:nt
                zprev, znow = CI[z1z2][1], CI[z1z2][2]
                flag = (isnan(R2[z1z2]) || R2[z1z2] < 0.99) ? "  ⚠" : ""
                @printf("  %2d→%-6d %11.6f %11.6f %9.4f %8d%s\n",
                    zprev, znow, Kfore_new[z1z2,1], Kfore_new[z1z2,2],
                    R2[z1z2], counts[z1z2], flag)
            end
            println("─"^62)
            @printf("Outer %2i | foredist = %.6f\n", outer_ct, foredist)
            @printf("K range for (%4.2f, %4.2f): %2.4f, %2.4f\n",
                policies.η, policies.τ,
                minimum(Kt[501:end]), maximum(Kt[501:end]))

            summarizeKtByTransition(Kt, it_t, params.π_z, burnin)
        end

        # damped update
        Kfore = λ_damp .* Kfore_new .+ (1 - λ_damp) .* Kfore
        outer_ct += 1
    end

    converged = foredist <= dTol
    if converged
        if verbose
            @printf("\nConverged in %i outer iters. foredist = %.6f\n", outer_ct-1, foredist)
        end
    else
        # always report failure, regardless of verbose
        @printf("\nHit maxout=%i without converging. foredist = %.6f\n", maxout, foredist)
        @printf("Maxout at policy (η, τ) = (%4.2f, %4.2f)\n", policies.η, policies.τ)
    end

    return Kfore, Kt, it_t
end


end # module