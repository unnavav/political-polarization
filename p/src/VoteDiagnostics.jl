# VoteDiagnostics.jl
# ALL OF THIS WAS WRITTEN BY CLAUDE, MINOR EDITS BY ME.
#
# Cross-sectional vs. aggregate decomposition of the vote differential
#
#   d(l,a; K,it) = V(K'_A; l,a) - V(K'_B; l,a)
#
# where K'_A, K'_B come from the two regime-specific KS forecast rules
# evaluated at the SAME (K, it). Splits d into
#
#   dbar(K,it) = Σ_{l,a} μ(l,a) d(l,a)          (common / aggregate part)
#   η(l,a)     = d(l,a) - dbar                   (household-specific part)
#
# The diagnostic ratio is  σ_η / sd(dbar),  where σ_η is the μ-weighted
# cross-sectional sd at a point and sd(dbar) is the variation of dbar
# across the aggregate states actually visited.
#
#   ratio >> 1  →  frozen: same households vote the same way forever
#   ratio ~ 1   →  live:   the marginal household flips over the cycle

module VoteDiagnostics

using Statistics: mean, std
using Printf: @printf, @sprintf 

using ..Compute: weight

export vote_diff, decompose_at, decompose_path, VoteDiagStats, log_diag


# ─── K-interpolation of V at a single (K, it) ───────────────────────────
#
# V is (nk, nt, nl, na); returns the (nl, na) slice at capital level Kval.
# Clamps to the grid ends and reports whether clamping occurred, since the
# regime-A rule gets evaluated at regime-B capital levels and vice versa.

function interpV(V::Array{Float64,4}, Kgrid::Vector{Float64},
                 Kval::Float64, it::Int)

    lo, hi = first(Kgrid), last(Kgrid)
    outside = (Kval < lo) || (Kval > hi)
    Kc = clamp(Kval, lo, hi)

    ix, we = weight(Kgrid, Kc)
    slice = we .* @view(V[ix, it, :, :]) .+ (1 - we) .* @view(V[ix+1, it, :, :])

    return slice, outside
end


# ─── forecast rule evaluation ───────────────────────────────────────────
#
# Kfore[it, :] = [const, coef on log K, (further coefs...)]
# Only the first two columns are used; extra columns are ignored so this
# still works if the rule later carries additional regressors.

function forecastK(Kfore::Matrix{Float64}, K::Float64, it::Int)
    return exp(Kfore[it, 1] + Kfore[it, 2] * log(K))
end


# ─── the vote differential itself ───────────────────────────────────────

"""
    vote_diff(V, Kgrid, KforeA, KforeB, K, it)

Returns `(d, KA, KB, n_outside)` where `d` is the (nl,na) array of
V(K'_A) - V(K'_B) and `n_outside` counts how many of the two projections
landed outside `Kgrid` (0, 1, or 2).
"""
function vote_diff(V::Array{Float64,4}, Kgrid::Vector{Float64},
                   KforeA::Matrix{Float64}, KforeB::Matrix{Float64},
                   K::Float64, it::Int)

    KA = forecastK(KforeA, K, it)
    KB = forecastK(KforeB, K, it)

    VA, outA = interpV(V, Kgrid, KA, it)
    VB, outB = interpV(V, Kgrid, KB, it)

    return VA .- VB, KA, KB, Int(outA) + Int(outB)
end


# ─── per-state decomposition ────────────────────────────────────────────

struct VoteDiagStats
    dbar::Float64        # μ-weighted mean of d
    sigma_eta::Float64   # μ-weighted cross-sectional sd of d
    theta::Float64       # μ-mass with d > 0
    near_mass::Float64   # μ-mass within ±band·σ_η of indifference
    KA::Float64
    KB::Float64
    n_outside::Int
end

"""
    decompose_at(V, Kgrid, KforeA, KforeB, μ, K, it; band=0.1)

Decomposition at a single aggregate state (K, it). `μ` is the (nl,na)
cross-sectional distribution; it is renormalized internally so it need not
sum to exactly 1.

`near_mass` is the share of households whose d lies within `band * σ_η` of
zero — the density-near-indifference term. A healthy σ_η/sd(dbar) ratio with
near-zero `near_mass` still gives a sluggish Θ, since dΘ/d(dbar) is governed
by that density.
"""
function decompose_at(V::Array{Float64,4}, Kgrid::Vector{Float64},
                      KforeA::Matrix{Float64}, KforeB::Matrix{Float64},
                      μ::Matrix{Float64}, K::Float64, it::Int;
                      band::Float64 = 0.1)

    d, KA, KB, nout = vote_diff(V, Kgrid, KforeA, KforeB, K, it)

    w = μ ./ sum(μ)

    dbar = sum(w .* d)
    var  = sum(w .* (d .- dbar).^2)
    σ_η  = sqrt(max(var, 0.0))

    θ = sum(w .* (d .> 0))

    near = σ_η > 0 ? sum(w .* (abs.(d) .<= band * σ_η)) : 0.0

    return VoteDiagStats(dbar, σ_η, θ, near, KA, KB, nout)
end


# ─── path aggregation ───────────────────────────────────────────────────

"""
    decompose_path(V, Kgrid, KforeA, KforeB, μs, Ks, its; band=0.1)

Runs `decompose_at` along a simulated path and assembles the summary the
diagnostic actually turns on.

`μs` may be either a Vector of (nl,na) distributions (one per period, the
honest version) or a single (nl,na) distribution reused at every point (the
cheap version, if you don't have the path of distributions stored).

Returns a NamedTuple with the per-period series plus:

  `ratio`      = mean(σ_η) / sd(dbar)   ← the headline number
  `sd_dbar`    = sd of dbar along the path (the "swing")
  `mean_sigma` = mean cross-sectional sd (the "spread")
"""
function decompose_path(V::Array{Float64,4}, Kgrid::Vector{Float64},
                        KforeA::Matrix{Float64}, KforeB::Matrix{Float64},
                        μs, Ks::Vector{Float64}, its::Vector{Int};
                        band::Float64 = 0.1)

    T = length(Ks)
    @assert length(its) == T "Ks and its must be the same length"

    getμ = μs isa AbstractVector ? (t -> μs[t]) : (t -> μs)

    dbar  = zeros(T); sigma = zeros(T); theta = zeros(T)
    near  = zeros(T); KAs   = zeros(T); KBs   = zeros(T)
    nout  = 0

    for t in 1:T
        s = decompose_at(V, Kgrid, KforeA, KforeB, getμ(t), Ks[t], its[t];
                         band = band)
        dbar[t]  = s.dbar
        sigma[t] = s.sigma_eta
        theta[t] = s.theta
        near[t]  = s.near_mass
        KAs[t]   = s.KA
        KBs[t]   = s.KB
        nout    += s.n_outside
    end

    sd_dbar    = T > 1 ? std(dbar) : 0.0
    mean_sigma = mean(sigma)
    ratio      = sd_dbar > 0 ? mean_sigma / sd_dbar : Inf

    return (dbar = dbar, sigma_eta = sigma, theta = theta, near_mass = near,
            KA = KAs, KB = KBs,
            sd_dbar = sd_dbar, mean_sigma = mean_sigma, ratio = ratio,
            frac_outside = nout / (2T),
            theta_min = minimum(theta), theta_max = maximum(theta),
            crossings = count(t -> (theta[t] - 0.5) * (theta[t+1] - 0.5) < 0,
                              1:T-1))
end


# ─── logging ────────────────────────────────────────────────────────────

"""
    log_diag(res; iter=nothing, label="")

One-line-per-block dump of `decompose_path` output, for calling once per
outer forecast iteration so the series can be watched as the rules update.
"""
function log_diag(res; iter = nothing, label::String = "")

    tag = iter === nothing ? label : @sprintf("iter %-4d %s", iter, label)

    @printf("\n─── vote decomposition %s ───\n", tag)
    @printf("  spread   mean σ_η      = %10.6f\n", res.mean_sigma)
    @printf("  swing    sd(dbar)      = %10.6f\n", res.sd_dbar)
    @printf("  RATIO    σ_η/sd(dbar)  = %10.4f   %s\n", res.ratio,
            res.ratio > 10 ? "(frozen)" :
            res.ratio > 3  ? "(sluggish)" : "(live)")
    @printf("  dbar     mean          = %10.6f\n", mean(res.dbar))
    @printf("  Θ        range         = [%6.4f, %6.4f]   crossings = %d\n",
            res.theta_min, res.theta_max, res.crossings)
    @printf("  mass near indifference = %10.6f\n", mean(res.near_mass))
    @printf("  K'_A     mean          = %10.4f\n", mean(res.KA))
    @printf("  K'_B     mean          = %10.4f\n", mean(res.KB))
    @printf("  frac projections off-grid = %7.4f\n", res.frac_outside)

    return nothing
end

end # module