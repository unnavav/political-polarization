module EGM

using ..ModelTypes: ModelParams, ImpliedRegimeParams, ProposedPolicies

using ..Compute: weight
using ..ModelFunctions: tax, u

# Numerical derivative of EV using spacing-weighted central differences.
# One-sided at boundaries. Clamped to [1e-12, 1e12].
function numdev(V0::Matrix{Float64}, agrid::Vector{Float64})
    nl, na = size(V0)
    DEV = zeros(nl, na)

    for ia in 1:na
        a = agrid[ia]

        if ia == 1
            DEV[:, ia] = (V0[:, ia+1] .- V0[:, ia]) ./ (agrid[ia+1] - a)
        elseif ia == na
            DEV[:, ia] = (V0[:, ia] .- V0[:, ia-1]) ./ (a - agrid[ia-1])
        else
            a0 = agrid[ia-1]
            a1 = agrid[ia+1]
            dr = (V0[:, ia+1] .- V0[:, ia]) ./ (a1 - a)
            dl = (V0[:, ia] .- V0[:, ia-1]) ./ (a - a0)
            # spacing-weighted average
            DEV[:, ia] = ((a - a0) / (a1 - a0)) .* dr .+
                         ((a1 - a) / (a1 - a0)) .* dl
        end
    end

    clamp!(DEV, 1e-12, 1e12)
    return DEV
end

# Single EGM iteration. Returns (TV, G, C).
function solve(V0::Array{Float64, 3}, terms::ModelParams, 
    policies::ProposedPolicies, prices::ImpliedRegimeParams, 
    CI::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}})

    lgrid  = terms.lgrid
    agrid  = terms.agrid
    zgrid  = terms.zgrid
    ϕ      = terms.ϕ
    σ      = terms.σ
    β      = terms.β
    α      = terms.α
    τ      = policies.τ
    λ      = prices.λ
    r      = prices.r
    w      = prices.w   
    captax = policies.captax

    nz = length(zgrid);

    nt, nl, na = size(V0)
    amin = agrid[1]
    amax = agrid[end]

    # --- Step 1: endogenous grid via Euler equation inversion ---
    D = zeros(nt, nl, na)
    for it = 1:nt
        Dt = numdev(V0[it, :, :], agrid);
        D[it,:, :] = Dt;
    end
    
    endoK = zeros(nt, nl, na)

    for it = 1:nt
        iz = CI[it][1];
        for ia in 1:na
            kpr = agrid[ia]
            for il in 1:nl
                c = (β * D[it, il, ia])^(-1.0 / σ)
                y = tax(w[iz] * lgrid[il], λ[iz], τ)
                numerator = c + kpr - y + r[iz] * (1.0 - captax[il]) * ϕ
                denom = 1.0 + r[iz] * (1.0 - captax[il])
                endoK[it, il, ia] = numerator / denom
            end
        end
    end

    # --- Step 2: interpolate decision rule onto exogenous grid ---
    G = zeros(nt, nl, na)

    for it = 1:nt
        for il in 1:nl
            lkvals = endoK[it, il, :]
            for ia in 1:na
                k = agrid[ia]
                if k <= lkvals[1]
                    G[it, il, ia] = amin
                elseif k > lkvals[end]
                    G[it, il, ia] = amax
                else
                    ix, we = weight(lkvals, k)
                    try
                        kpr = we * agrid[ix] + (1.0 - we) * agrid[ix+1]
                        G[it, il, ia] = max(amin, kpr)
                    catch
                        error("Interpolation error in EGM solve.")
                        display("it = $it, il = $il, ia = $ia, k = $k, ix = $ix, we = $we")
                    end
                end
            end
        end
    end

    clamp!(G, amin, amax)

    # --- Step 3: consumption and value function update ---
    C  = zeros(nt, nl, na)
    TV = zeros(nt, nl, na)

    for it = 1:nt
        iz = CI[it][1];
        for ia in 1:na
            for il in 1:nl
                c = (1.0 + r[iz] * (1.0 - captax[il])) * agrid[ia] +
                    tax(w[iz] * lgrid[il], λ[iz], τ) - G[it, il, ia] -
                    r[iz] * (1.0 - captax[il]) * ϕ
                C[it, il, ia] = max(1e-6, c)

                ix, we = weight(agrid, G[it, il, ia])
                ev = we * V0[it, il, ix] + (1.0 - we) * V0[it, il, ix+1]
                TV[it, il, ia] = u(C[it, il, ia], σ) + β * ev
            end
        end
    end

    return TV, G, C
end

# One-sided finite difference for a 1D expected value vector.
function solveD(EV0::Vector{Float64}, ik::Int, kchgrid::Vector{Float64})
    nk = length(kchgrid)
    k = kchgrid[ik]
    if ik == nk
        return (EV0[ik] - EV0[ik-1]) / (k - kchgrid[ik-1])
    else
        return (EV0[ik+1] - EV0[ik]) / (kchgrid[ik+1] - k)
    end
end

end # module