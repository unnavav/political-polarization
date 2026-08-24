# ModelFunctions.jl
module ModelFunctions

using ..Compute: weight, weight_vec
using ..ModelTypes: ModelParams, ImpliedRegimeParams, ProposedPolicies
using ..DistrTools: getDistr, transitDistr
using LinearAlgebra: dot
export u, bellmanValue, getExpectationKS, tax, mapVotes, calcr, calcw, dwde, drde, mapVotesKS, voteSharePath, build_votes


# ─── Utility fns ───

# CRRA utility
function u(c::Float64, σ::Float64)
    c = max(c, 1e-12)
    if σ == 1.0
        return log(c)
    else
        return c^(1.0 - σ) / (1.0 - σ)
    end
end

# Bellman utility (consumption from budget + continuation)
function bellmanValue(apr::Float64, y::Float64, β::Float64,
              σ::Float64, vc::Float64)
    if σ == 1.0
        return log(y - apr) + β * vc
    else
        return (y - apr)^(1.0 - σ) / (1.0 - σ) + β * vc
    end
end

function getExpectationKS(futureKs::Matrix{Float64}, V0::Array{Float64,4}, params::ModelParams,
                    CI::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}},
                    LI::LinearIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}})
    nk, nt, nl, na = size(V0)
    π_z = params.π_z
    π_l = params.π_l
    Kgrid = params.Kgrid
    EV = zeros(nk, nt, nl, na)
    
    nz = length(params.zgrid);

    @views for ik in 1:nk
        for it in 1:nt
            EK = futureKs[it, ik]
            ix, we = weight(Kgrid, EK)
            iz = CI[it][2]
            for il in 1:nl
                for ia in 1:na
                    ev = 0.0
                    for jz in 1:nz
                        jt = LI[iz, jz];
                        weighted_V = we*V0[ix, jt, :, ia] + (1-we)*V0[ix+1, jt, :, ia]
                        ev += π_z[iz, jz]*dot(π_l[il, :], weighted_V)
                    end
                    EV[ik, it, il, ia] = ev
                end
            end
        end
    end

    return EV
end

# ─── Government and Elections ───

# HSV tax function
function tax(gross::Float64, λ::Float64, τ::Float64)
    return λ * gross^(1.0 - τ)
end


function mapVotes(VOTES::Array{Float64,2}, params, μ_today::Array{Float64,2})

    agrid = params.agrid; amu = params.amu;

    nl, _ = size(VOTES); nμ = length(amu); 

    ixgrid = zeros(Int, nμ)
    wegrid = zeros(Float64, nμ)
    for iμ in 1:nμ
        ix, we = weight(agrid, amu[iμ])
        ixgrid[iμ] = ix
        wegrid[iμ] = we
    end

    # interpolation sanity: each row's weights must sum to 1 (mass-preserving)
    for im in 1:nμ
        @assert isapprox(wegrid[im] + (1 - wegrid[im]), 1.0; atol=1e-10) "weights don't sum to 1 at im=$im"
    end

    vdistr = zeros(nl, nμ)

    for il in 1:nl, iμ in 1:nμ
        μ_val = μ_today[il, iμ]
        if μ_val > 0.0
            ix = ixgrid[iμ]
            we = wegrid[iμ]
            vdistr[il, iμ] += VOTES[il, ix] * we * μ_val +
                                VOTES[il, ix+1] *(1- we) * μ_val
        end
    end

    voteshare = sum(vdistr)

    @assert -1e-8 ≤ voteshare ≤ 1.0 + 1e-8 "Θ out of [0,1]: $voteshare"

    winner = voteshare > 0.5 ? 1 : 0

    return vdistr, voteshare 

end

function mapVotesKS(votes::Array{Int,3}, μ::Array{Float64,3},
                    futureK::Float64, Kgrid::Vector{Float64},
                    amu::Vector{Float64}, agrid::Vector{Float64})
    nl, nmu = size(μ, 2), size(μ, 3)
    ixK, weK = weight(Kgrid, futureK)
    todays_Vote = weK .* votes[ixK, :, :] .+ (1-weK) .* votes[ixK+1, :, :]  # (nl,na)

    vdistr = zeros(nl, nmu)
    share = 0.0
    for im in 1:nmu
        ix, we = weight(agrid, amu[im])
        for il in 1:nl
            mass = μ[1, il, im]
            v = we*todays_Vote[il, ix] + (1-we)*todays_Vote[il, ix+1]
            vdistr[il, im] = mass * v
            share += mass * v
        end
    end

    return vdistr, share
end

function getVotes(V0::Array{Float64,5}, params::ModelParams, futureKs::Matrix{Float64}, CI, LI)
    nθ, nk, nt, nl, na = size(V0)
    EV = zeros(Float64, nθ, nk, nt, nl, na);

    # getting transition probabilities etc
    kernel = params.kernel; 
    π_Θ = params.π_Θ; #this has many jobs, but for now let's just call it the transition matrix
    ε_mean = mean(kernel); 
    if abs(ε_mean) > 1e-8
        @warn "Kernel mean is not ~0, but $ε_mean. Check the kernel distribution and residuals."
    end


    for iθ in 1:nθ
        EV[iθ, :, :, :, :] = getExpectationKS(futureKs, V0[iθ, :, :, :, :], params, CI, LI);
    end

    # here's the problem. Our nθ is coarse, kernel is not. So I'll need to interpolate EV to be
    # more granular.  But since this is, at its heart, a discrete problem, I can't just do a cubic
    # spline across all of EV. I've logspaced around 0.5, so what I'll need to do is interpolate the 
    # Θgrid on the kernel grid, and then use that to get an EV that matches our kernel. 

    # CLAUDE: IF i COPY-PASTE THIS, DOUBLE CHECK THAT THE KERNEL IS STILL LINEARLY SPACED
    n_knots = (length(kernel)-1)÷2+1;
    knots = LinRange(0, 1, n_knots)

    ixs, wes = weight_vec(params.Θgrid, collect(knots));
    EV_k = zeros(n_knots, nk, nt, nl, na);
    for ik in 1:n_knots
        ix = ixs[ik]; we = wes[ik];
        EV_k[ik, :, :, :, :] = we * EV[ix, :, :, :, :] + (1-we)*EV[ix+1, :, :, :, :];
    end

    # finally computing votes 😭🫠🚬
    votes = zeros(Int, n_knots, nk, nt, nl, na);
    ix, we = weight(knots, 0.5); #we should be 1 here; the knot should be already 0.5
    @assert isapprox(we, 1.0; atol=1e-10) "0.5 is not exactly a knot: n_knots=$n_knots gives we=$we. Need odd n_knots."
    A_range = ix:n_knots    # Θ' >= 0.5 → Incumbent
    B_range = 1:ix-1            # Θ' < 0.5 → Challenger
    PROB_TOL = 1e-10      # just making sure this doesn't return NA or blow up

    for ik in 1:n_knots
        a_prob = sum(π_Θ[ik, A_range])
        b_prob = sum(π_Θ[ik, B_range])

        a_reachable = a_prob > PROB_TOL
        b_reachable = b_prob > PROB_TOL

        for iK in 1:nk, it in 1:nt, il in 1:nl, ia in 1:na
            if a_reachable && b_reachable
                W_A = dot(π_Θ[ik, A_range], @view EV_k[A_range, iK, it, il, ia]) / a_prob
                W_B = dot(π_Θ[ik, B_range], @view EV_k[B_range, iK, it, il, ia]) / b_prob
                votes[ik, iK, it, il, ia] = W_A > W_B
            elseif a_reachable
                votes[ik, iK, it, il, ia] = 1      # only A reachable → vote A
            else
                votes[ik, iK, it, il, ia] = 0      # only B reachable (or neither) → vote B
            end
        end
    end

    return votes
end

function voteSharePath(incumbent, challenger, zt_shock, params)
    Kgrid = params.Kgrid
    amu   = params.amu
    agrid = params.agrid
    nk, nz, nl, na = size(incumbent.V)
    NT = length(zt_shock)

    # --- 1. Build the vote indicator votes[ik, il, ia] once ---
    # EV for each policy (full, not collapsed — we need every K-node)
    fK_i = exp.(incumbent.Kfore[:,1]  .+ incumbent.Kfore[:,2]  .* log.(Kgrid)')
    fK_c = exp.(challenger.Kfore[:,1] .+ challenger.Kfore[:,2] .* log.(Kgrid)')
    EV_i = getExpectationKS(fK_i, incumbent.V,  params)   # (nk,nz,nl,na)
    EV_c = getExpectationKS(fK_c, challenger.V, params)

    # --- 2. Simulate the distribution forward along zt_shock ---
    # seed μ from a starting K (grid median) and z = zt_shock[1]
    ik0 = cld(nk, 2)
    Kt  = zeros(NT+1); Kt[1] = Kgrid[ik0]

    # initial cross-section: use the incumbent's policy at the seed
    # (reuse your getDistr → collapse → slice machinery, or seed uniform and burn in)
	CI_dist = CartesianIndices(params.π_z)   # (nz,nz) for getDistr
	LI_dist = LinearIndices(params.π_z)
    G_start = incumbent.G[ik0, :, :, :]                        # (nz,nl,na) single-z
	G_pair = zeros(nz*nz, nl, na)
	for it in 1:(nz*nz)
		today = CI_dist[it][2]
		G_pair[it, :, :] = G_start[today, :, :]
	end

    # ... lift to pair-state (nz²) if getDistr is pair-state ...
    μ_transit, _ = getDistr(G_pair, params.amu, params.agrid, params.π_l, params.π_z,
                         CI_dist, LI_dist, params.ϕ)

    nmu = length(params.amu)
    μ_today = zeros(nz, nl, nmu)
    for it in 1:(nz*nz)
        today = CI_dist[it][2]
        μ_today[today, :, :] .+= μ_transit[it, :, :]
    end
    μ_transit = μ_today[zt_shock[1]:zt_shock[1], :, :]   # (1, nl, nmu)

    share_path = zeros(NT)

    for t in 1:NT
        K  = Kt[t]
        iz = zt_shock[t]

        # forecasted future K from the INCUMBENT rule at today's (K, z)
        futureK = exp(incumbent.Kfore[iz,1] + incumbent.Kfore[iz,2]*log(K))

        # vote share at this state: z-slice the indicator, interp to futureK
        votes_z = build_votes(EV_c, EV_i, iz)            # (nl, na) → see below
        _, share_path[t] = mapVotesKS(votes_z, μ_transit, futureK, Kgrid, amu, agrid)

        # step the distribution forward with the incumbent's policy at (K, z)
        ix, we = weight(Kgrid, K)
        G_t = we .* incumbent.G[ix, iz, :, :] .+ (1-we) .* incumbent.G[ix+1, iz, :, :]
        μ_transit, Kt[t+1] = transitDistr(G_t, μ_transit, amu, agrid, params.ϕ, params.π_l)
    end

    return share_path, Kt
end

function build_votes(EV_c::Array{Float64,4}, EV_i::Array{Float64,4}, it::Int)
    return Int.(EV_c[:, it, :, :] .> EV_i[:, it, :, :])
end

function build_transition_matrix(kernel, ρ0::Float64, ρ1::Float64)

    n_knots = (length(kernel)-1)÷2+1;
    knots = LinRange(0, 1, n_knots);

    π_Θ = zeros(n_knots, n_knots);
    for i in 1:n_knots
        m = ρ0 + ρ1 * knots[i]              # where I guess I'll be tomorrow given this θ
        ix, we = weight(knots, m)                    # m sits between knot ix and ix+1
        
        # this is insane but let me cook i swear

        # rolling window of length n_knots into the length-(2n_knots-1) kernel.
        # ix knots below, (n_knots - ix) above, for the lower placement:
        lo_start = n_knots - ix + 1             # so the window is exactly n_knots long
        lo = lo_start : lo_start + n_knots - 1
        hi = min.(lo .+ 1, length(kernel))     # the adjacent placement, shifted by one
    
        π_Θ[i, :] = we * kernel[lo] .+ (1-we) * kernel[hi];

        s = sum(π_Θ[i, :])
        if s > 0
            π_Θ[i, :] ./= s                      # renormalize after truncation
        end
    end

    return π_Θ
end

# ─── Pricing and Policy Responsiveness ───

# Competitive prices from firm FOCs (with migration rate η)
function calcr(α::Float64, δ::Float64, k::Float64, η::Float64,z::Vector{Float64})
    return α .* z .* (k / (1.0 + η)).^(α - 1.0) .- δ
end

function calcr(α, δ, k::Float64, η, z::Float64) #scalar method
     return α * z * (k/(1+η))^(α-1) - δ
end

function calcw(α::Float64, k::Float64, η::Float64,z::Vector{Float64})
    return (1.0 - α) .* z .* (k / (1.0 + η)).^α
end

function calcw(α, k::Float64, η, z::Float64) #scalar method
    return (1-α) * z * (k/(1+η))^α
end

# Partials w.r.t. η (for comparative statics / GE effects)
function dwde(α::Float64, k::Float64, η::Float64)
    return -α * (1.0 - α) * k^α * (1.0 + η)^(-α - 1.0)
end

function drde(α::Float64, k::Float64, η::Float64)
    return α * (1.0 - α) * k^(α - 1.0) * (1.0 + η)^(-α - 1.0)
end

end