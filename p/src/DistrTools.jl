module DistrTools

export getDistr, transitDistr, condense

using ..Compute: weight, supnorm
using LinearAlgebra: dot
using Printf: @printf

# ─── Stationary distribution ───

function getDistr(G::Array{Float64,3}, amu::Vector{Float64}, agrid::Vector{Float64},
                 π_l::Matrix{Float64}, π_z::Matrix{Float64},
                 CI::CartesianIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}},
                 LI::LinearIndices{2, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}},
                 ϕ::Float64; verbose::Bool=false,
                 vTol::Float64=1e-6)

    nt, nl, _ = size(G)
    nmu = length(amu)
    nz = size(π_z, 1)

    ixgrid::Array{Int,3} = zeros(nt, nl, nmu)
    wegrid::Array{Float64,3} = zeros(nt, nl, nmu)

    # precompute interpolation maps
    for it in 1:nt
        for im in 1:nmu
            kval = amu[im]
            for il in 1:nl
                ix, we = weight(agrid, kval)
                kdval = G[it, il, ix] * we + G[it, il, ix+1] * (1.0 - we)
                kdval = clamp(kdval, amu[1], amu[end])       # ← add this line

                ix2, we2 = weight(amu, kdval)
                ixgrid[it, il, im] = ix2
                wegrid[it, il, im] = we2
            end
        end
    end

    μ = ones(nt, nl, nmu) / (nt * nmu * nl)
    distance = 20.0
    iter_ct = 1

    while distance > vTol
        μ1 = zeros(nt, nl, nmu)

        for it in 1:nt, im in 1:nmu, il in 1:nl
            ix = ixgrid[it, il, im]
            we = wegrid[it, il, im]
            μ_val = μ[it, il, im]
            today = CI[it][2]
            if μ_val > 0.0
                for tomorrow in 1:nz
                    jt = LI[today, tomorrow]
                    for jl in 1:nl
                        base = π_z[today, tomorrow] * π_l[il, jl] * μ_val
                        if ix < nmu
                            μ1[jt, jl, ix]     += base * we
                            μ1[jt, jl, ix + 1] += base * (1.0 - we)
                        else
                            μ1[jt, jl, ix]     += base          # all mass at last node, nothing dropped
                        end
                    end
                end
            end
        end

        distance = supnorm(μ1, μ, 3)
        iter_ct += 1
        μ .= μ1
    end

    if verbose
        @printf("\n\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f\n",
                iter_ct, distance, sum(μ))
    end

    distrAgg = dropdims(sum(μ, dims=(1, 2)), dims=(1, 2))   # (nmu,)
    kagg = dot(amu .- ϕ, distrAgg)

    return μ, kagg
end

# ─── Transition distribution (one-period forward) ───

function transitDistr(g_t::Matrix{Float64}, μ_prev::Array{Float64,3},
                      amu::Vector{Float64}, agrid::Vector{Float64},
                      ϕ::Float64, pil::Matrix{Float64})

    nl, _ = size(g_t)
    nd, _, nmu = size(μ_prev)

    ixgrid = zeros(Int, nl, nmu)
    wegrid = zeros(nl, nmu)

    for im in 1:nmu
        kval = amu[im]
        for il in 1:nl
            ix, we = weight(agrid, kval)
            kdval = g_t[il, ix] * we + g_t[il, ix+1] * (1.0 - we)
            kdval = clamp(kdval, amu[1], amu[end])
            ix2, we2 = weight(amu, kdval)
            ixgrid[il, im] = ix2
            wegrid[il, im] = we2
        end
    end

    μ1 = zeros(size(μ_prev))

    for id in 1:nd, im in 1:nmu, il in 1:nl
        ix = ixgrid[il, im]
        we = wegrid[il, im]
        muval = μ_prev[id, il, im]

        if muval > 0.0
            for jl in 1:nl
                base = pil[il, jl] * muval
                if ix < nmu
                    μ1[id, jl, ix]     += base * we
                    μ1[id, jl, ix + 1] += base * (1.0 - we)
                else
                    μ1[id, jl, ix]     += base          # all mass at last node
                end
            end
        end
    end

    distrAgg = dropdims(sum(μ1, dims=(1, 2)), dims=(1, 2))
    kagg = dot(amu .- ϕ, distrAgg)

    return μ1, kagg
end

# ─── Distribution condensing ───

function condense(adistr::Array{Float64,3}, amu::Vector{Float64}, agrid::Vector{Float64})
    nd, nl, nmu = size(adistr)
    na = length(agrid)
    distr = zeros(nd, nl, na)

    for id in 1:nd, im in 1:nmu
        ix, we = weight(agrid, amu[im])
        distr[id, :, ix]   .+= we .* adistr[id, :, im]
        distr[id, :, ix+1] .+= (1.0 - we) .* adistr[id, :, im]
    end

    return distr
end

end