# this file contains all the main computational functions for the model, including interpolation, grid construction, and VFI with GSS. 
# Note that this file does not depend on any of the other files in the project, except for ModelTypes.jl, which defines the ModelParams
# struct. This allows us to avoid circular dependencies and keep the code organized.
# may 2026
# vaasavi, translated to Julia by Claude

module Compute

using Distributions: Normal, cdf
using Printf
using Statistics: mean, std, median

export weight, linterpolate, gss, getKgrid, getTauchen, dist, logspace, stationary, supnorm, summarizeKtByTransition

# ─── Interpolation ───

function weight(grid::Vector{Float64}, f::Float64)
    n = length(grid)
    if f >= grid[end]
        ix = n - 1
        we = 0.0
    elseif f < grid[1]
        ix = 1
        we = 1.0
    else
        ix = searchsortedlast(grid, f)
        we = (grid[ix+1] - f) / (grid[ix+1] - grid[ix])
    end
    return ix, we
end

function linterpolate(Vvec::Vector{Float64}, grid::Vector{Float64}, vi::Float64)
    n = length(grid)
    if vi <= grid[1]
        return Vvec[1]
    elseif vi >= grid[end]
        return Vvec[n]
    else
        il = searchsortedlast(grid, vi)
        wl = (grid[il+1] - vi) / (grid[il+1] - grid[il])
        return wl * Vvec[il] + (1.0 - wl) * Vvec[il+1]
    end
end

# ─── Grid construction ───

function getKgrid(nk::Int, kl::Float64, kh::Float64)
    raw = 10.0 .^ range(log10(1.0), log10(kh - kl + 1.0), length=nk)
    return raw .+ (kl - 1.0)
end

function logspace(l::Float64, h::Float64, n::Int)
    grid = exp10.(range(log10(1), log10(h-l+1), length=n));
    return grid .+ (l - 1.0)
end

# ─── Tauchen discretization ───

function getTauchen(Nz::Int, mu::Float64, sigma::Float64, rho::Float64, s::Float64)

    # need to back out σ^2_e given σ^2_l
    sigma_x = sqrt(sigma^2 / (1.0 - rho^2))

    x_1  = mu - s * sigma_x
    x_Nz = mu + s * sigma_x
    x_grid = range(x_1, x_Nz, length=Nz) |> collect
    z_grid = exp.(x_grid)

    w = x_grid[2] - x_grid[1]
    nd = Normal(mu, sigma)

    P_mat = zeros(Nz, Nz)

    for r in 1:Nz
        x_curr = x_grid[r] * rho

        P_mat[r, 1]  = cdf(nd, x_grid[1]  - x_curr + w / 2.0)
        P_mat[r, Nz] = 1.0 - cdf(nd, x_grid[Nz] - x_curr - w / 2.0)

        for c in 2:(Nz-1)
            upper = cdf(nd, x_grid[c] - x_curr + w / 2.0)
            lower = cdf(nd, x_grid[c] - x_curr - w / 2.0)
            P_mat[r, c] = upper - lower
        end
    end

    return P_mat, z_grid
end

function stationary(P::Matrix{Float64}; tol=1e-12, maxiter=10000)
    n = size(P, 1)
    π = ones(n) / n              # uniform start

    for _ in 1:maxiter
        π_new = P' * π           # one step forward
        if maximum(abs.(π_new - π)) < tol
            return π_new
        end
        π = π_new
    end
    error("didn't converge")
end

# ─── Sup-norm distance ───

function supnorm(M::Array{Float64}, N::Array{Float64}, nd::Int)
    x = abs.(M .- N)
    for _ in 1:nd
        x = maximum(x, dims=1)
    end
    return x[1]
end

# ─── Summary Statistics ───

function summarizeKtByTransition(Kt, it_t, π_z, burn_in)
    CI = CartesianIndices(π_z)
    nt = length(LinearIndices(π_z))

    println("\nKt summary by transition (post burn-in)")
    println("─"^72)
    @printf("  %-8s %8s %10s %10s %10s %10s %10s\n",
            "z₋₁→z", "n", "mean", "std", "min", "median", "max")
    println("─"^72)

    for it in 1:nt
        zprev, znow = CI[it][1], CI[it][2]
        # periods (post burn-in, with valid t+1) whose pair == it
        idx = [t for t in (burn_in+1):(length(Kt)-1) if it_t[t] == it]

        if isempty(idx)
            @printf("  %2d→%-5d %8d %10s %10s %10s %10s %10s\n",
                    zprev, znow, 0, "—", "—", "—", "—", "—")
            continue
        end

        vals = Kt[idx]
        @printf("  %2d→%-5d %8d %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                zprev, znow, length(idx),
                mean(vals), std(vals), minimum(vals), median(vals), maximum(vals))
    end
    println("─"^72)

    # overall, for reference
    allidx = (burn_in+1):(length(Kt)-1)
    allvals = Kt[allidx]
    @printf("  %-8s %8d %10.4f %10.4f %10.4f %10.4f %10.4f\n",
            "ALL", length(allvals), mean(allvals), std(allvals),
            minimum(allvals), median(allvals), maximum(allvals))
    println("─"^72)
end



end  # module