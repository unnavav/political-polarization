# ModelFunctions.jl
module ModelFunctions

using ..Compute: weight
export u, bellmanValue, tax, mapVotes, calcr, calcw, dwde, drde


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

# ─── Government and Elections ───

# HSV tax function
function tax(gross::Float64, λ::Float64, τ::Float64)
    return λ * gross^(1.0 - τ)
end


function mapVotes(VOTES::Array{Float64,3}, amu::Vector{Float64},
                   agrid::Vector{Float64}, adistr::Array{Float64,3},
                   pctDem::Float64)

    nl, np, nm = size(adistr)

    ixgrid = zeros(Int, nm)
    wegrid = zeros(nm)
    for im in 1:nm
        ix, we = weight(agrid, amu[im])
        ixgrid[im] = ix
        wegrid[im] = we
    end

    vdistr = zeros(nl, nm)

    for im in 1:nm, ip in 1:np, il in 1:nl
        muval = adistr[il, ip, im]
        if muval > 0.0
            ix = ixgrid[im]
            we = wegrid[im]
            vdistr[il, im] += pctDem     * VOTES[il, ix,   1] * we     * muval +
                              (1-pctDem) * VOTES[il, ix,   2] * we     * muval +
                              pctDem     * VOTES[il, ix+1, 1] * (1-we) * muval +
                              (1-pctDem) * VOTES[il, ix+1, 2] * (1-we) * muval
        end
    end

    majority = sum(vdistr)
    winner = majority > 0.5 ? 1 : 0

    return vdistr, winner
end

# ─── Pricing and Policy Responsiveness ───

# Competitive prices from firm FOCs (with migration rate η)
function calcr(α::Float64, δ::Float64, k::Float64, η::Float64,z::Vector{Float64})
    return α .* z .* (k / (1.0 + η)).^(α - 1.0) .- δ
end

function calcw(α::Float64, k::Float64, η::Float64,z::Vector{Float64})
    return (1.0 - α) .* z .* (k / (1.0 + η)).^α
end

# Partials w.r.t. η (for comparative statics / GE effects)
function dwde(α::Float64, k::Float64, η::Float64)
    return -α * (1.0 - α) * k^α * (1.0 + η)^(-α - 1.0)
end

function drde(α::Float64, k::Float64, η::Float64)
    return α * (1.0 - α) * k^(α - 1.0) * (1.0 + η)^(-α - 1.0)
end

end