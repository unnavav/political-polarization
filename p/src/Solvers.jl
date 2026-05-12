module Solvers

using LinearAlgebra: Tridiagonal

using ..ModelTypes: ModelParams, ImpliedRegimeParams, ProposedPolicies
using ..Compute: weight, supnorm
using ..EGM: solve
using ..ModelFunctions: tax, u, bellmanValue

using Printf: @printf

# ─── Golden section search ───

function gss(Cvals::Vector{Float64}, y::Float64, beta::Float64,
             sigma::Float64, searchgrid::Vector{Float64}, prec::Float64)

    r = (3.0 - sqrt(5.0)) / 2.0

    a = searchgrid[1]
    b = min(searchgrid[end], y)
    c = (1.0 - r) * a + r * b
    d = r * a + (1.0 - r) * b

    vc = linterpolate(Cvals, searchgrid, c)
    fc = -bellmanValue(c, y, beta, sigma, vc)

    vd = linterpolate(Cvals, searchgrid, d)
    fd = -bellmanValue(d, y, beta, sigma, vd)

    while abs(a - b) > prec
        if fc > fd
            a = c
            c = d
            d = r * a + (1.0 - r) * b
            fc = fd
            vd = linterpolate(Cvals, searchgrid, d)
            fd = -bellmanValue(d, y, beta, sigma, vd)
        else
            b = d
            d = c
            c = (1.0 - r) * a + r * b
            fd = fc
            vc = linterpolate(Cvals, searchgrid, c)
            fc = -bellmanValue(c, y, beta, sigma, vc)
        end
    end

    aval = c
    vc = linterpolate(Cvals, searchgrid, c)
    res = bellmanValue(c, y, beta, sigma, vc)
    return res, aval
end


# ─── VFI solver (calls EGM internally in MATLAB; loop structure preserved) ───

function EGMsolve(nl::Int, na::Int, terms::ModelParams, policies::ProposedPolicies, vTol::Float64; verbose::Bool=false)

    lgrid  = terms.lgrid
    agrid  = terms.agrid
    pil    = terms.pil
    
    r = terms.r
    w = terms.w

    τ = policies.τ
    λ = policies.λ
    η = policies.η
    captax = policies.captax

    V  = zeros(nl, na)
    G  = zeros(nl, na)
    TG = zeros(nl, na)
    TV = zeros(nl, na)
    V0 = zeros(nl, na)

    # warm-start V with log income
    scale = 0.25
    for ia in 1:na
        kval = agrid[ia]
        for il in 1:nl
            yval = scale * (1.0 + r * (1.0 - captax[il])) * kval +
                   w * lgrid[il] - r * phi
            V[il, ia] = log(max(1e-10, yval))
        end
    end

    # initial expected values
    for ia in 1:na
        for il in 1:nl
            V0[il, ia] = dot(pil[il, :], V[:, ia])
        end
    end

    d = 1e5
    iter_ct = 1

    while d > vTol
        # EGM step happens here — TV, G, C are set by your EGM.solve call
        # (In MATLAB this was implicit via the loop; you'll wire EGM.solve in)

        # here is where I need to establish the statespace first, and come back and fix the EGM. 

        d    = supnorm(TV, V, 2)
        kdst = supnorm(TG, G, 2)

        iter_ct += 1
        V  .= TV
        TG .= G

        for ia in 1:na
            for il in 1:nl
                V0[il, ia] = dot(pil[il, :], V[:, ia])
            end
        end
    end

    if verbose
        @printf("\n\tIteration %i: ||TV - V|| = %4.6f\t||TG - G|| = %4.6f\n",
                iter_ct, d, kdst)
    end

    return V, G, similar(G), V0   # (V, G, C, V0)
end

# ─── Cubic spline ───

function cubicSpline(l::Float64, h::Float64, r::Int, VK::Vector{Float64})
    val_grid = 10.0 .^ range(log10(l), log10(h), length=r+2) |> collect

    fti = zeros(r + 2)
    dti = zeros(r + 2)
    ttf = zeros(r + 2)

    fti[1] = VK[1]

    for i in 2:(r+2)
        dti[i] = val_grid[i] - val_grid[i-1]
        fti[i] = VK[i]
        ttf[i] = (fti[i] - fti[i-1]) / dti[i]
    end

    upper_diag = zeros(r + 1)
    lower_diag = zeros(r + 1)
    principal  = zeros(r + 1)

    for i in 2:r
        upper_diag[i] = dti[i]
        lower_diag[i] = dti[i+2]
    end

    upper_diag[2] += dti[2]^2 / dti[3]
    lower_diag[r] += dti[r]^2 / dti[r-1]

    omega_1 = dti[3] - dti[2]^2 / dti[3]
    omega_r = dti[r+1] - dti[r+2]^2 / dti[r+1]

    principal[1]   = 2.0 * (dti[2] + dti[3]) - omega_1
    principal[r+1] = 2.0 * (dti[r+1] + dti[r+2]) - omega_r

    for i in 2:r
        principal[i] = 2.0 * (dti[i+1] + dti[i+2])
    end

    # RHS
    f = zeros(r)
    f[1] = 3.0 * (dti[3] * ttf[2] + dti[2] * ttf[3]) -
           2.0 * (dti[3] * ttf[2] - dti[2]^2 / dti[3] * ttf[3])
    f[r] = 3.0 * (dti[r+2] * ttf[r+1] + dti[r+1] * ttf[r+2]) -
           2.0 * (dti[r+1] * ttf[r+2] - dti[r+2]^2 / dti[r+1] * ttf[r+1])

    for i in 2:r
        f[i] = 3.0 * (dti[i+2] * ttf[i+1] + dti[i+1] * ttf[i+2])
    end

    # Build tridiagonal and solve
    dl = [lower_diag[i] for i in 2:r]      # sub-diagonal, length r-1
    dp = [principal[i]  for i in 1:r]
    du = [upper_diag[i] for i in 2:r]       # super-diagonal, length r-1
    # Adjust lengths for Tridiagonal constructor
    T = Tridiagonal(dl[1:r-1], dp, du[1:r-1])
    s = T \ f

    s_0 = 2.0 * ttf[2] - (dti[2] / dti[3])^2 * ttf[3] -
          (1.0 - (dti[2] / dti[3])^2) * s[1] +
          (dti[2] / dti[3])^2 * s[2]

    s_r1 = 2.0 * (ttf[r+2] - (dti[r+2] / dti[r+1])^2 * ttf[r+1]) -
           (1.0 - (dti[r+2] / dti[r+1])^2) * s[r] +
           (dti[r+2] / dti[r+1])^2 * s[r-1]

    s_fin = [s_0; s; s_r1]

    C = zeros(r + 1, 4)
    for i in 1:(r+1)
        C[i, 1] = fti[i]
        C[i, 2] = s_fin[i]
        C[i, 3] = 3.0 * ttf[i+1] / dti[i+1] -
                   2.0 * s_fin[i] / dti[i+1] - s_fin[i+1] / dti[i+1]
        C[i, 4] = (-2.0 * ttf[i+1] + s_fin[i] + s_fin[i+1]) / dti[i+1]^2
    end

    return C
end

function getSplineVal(coeffs::Matrix{Float64}, lval::Float64, lgrid::Vector{Float64})
    il = searchsortedlast(lgrid, lval)
    if il == length(lgrid)
        il -= 1
    end
    d = lval - lgrid[il]
    return coeffs[il, 1] + coeffs[il, 2] * d +
           coeffs[il, 3] * d^2 + coeffs[il, 4] * d^3
end

# ─── VFI with GSS (interpV and backsolve) ───

function interpsolve(terms::NamedTuple, V::Matrix{Float64},
                 EV::Matrix{Float64}, vTol::Float64)

    agrid  = terms.agrid;  lgrid = terms.lgrid
    na = length(agrid);    nl = length(lgrid)
    beta = terms.beta;     sigma = terms.sigma
    r = terms.r;           w = terms.w
    lambda = terms.lamval; tau = terms.tau
    captax = terms.captax; pil = terms.pil

    TV = zeros(nl, na)
    G  = zeros(nl, na)

    for il in 1:nl
        l = lgrid[il]
        tau_r = captax[il]
        for ia in 1:na
            a = agrid[ia]
            y = tax(w * l, lambda, tau) + (1.0 + r * (1.0 - tau_r)) * a
            ix = searchsortedlast(agrid, y)
            ix = max(ix, 1)
            searchgrid = agrid[1:ix]
            vval, aval = gss(EV[il, :], y, beta, sigma, searchgrid, vTol * 1e-2)
            TV[il, ia] = vval
            G[il, ia]  = aval
        end
    end

    for ia in 1:na, il in 1:nl
        EV[il, ia] = dot(pil[il, :], TV[:, ia])
    end

    return TV, G, EV
end

function backsolve(terms::NamedTuple, EV::Matrix{Float64}, vTol::Float64)

    agrid  = terms.agrid;  lgrid = terms.lgrid
    na = length(agrid);    nl = length(lgrid)
    beta = terms.beta;     sigma = terms.sigma
    r = terms.r;           w = terms.w
    lambda = terms.lamval; tau = terms.tau
    captax = terms.captax

    V = zeros(nl, na)
    G = zeros(nl, na)

    for il in 1:nl
        l = lgrid[il]
        tau_r = captax[il]
        for ia in 1:na
            a = agrid[ia]
            y = tax(w * l, lambda, tau) + (1.0 + r * (1.0 - tau_r)) * a
            ix = searchsortedlast(agrid, y)
            ix = max(ix, 1)
            searchgrid = agrid[1:ix]
            vval, aval = gss(EV[il, :], y, beta, sigma, searchgrid, vTol * 1e-2)
            V[il, ia] = vval
            G[il, ia] = aval
        end
    end

    return V, G
end

# ─── Backsolve: GSS-based solve given a fixed continuation value ───

function GSSbacksolve(nl::Int, na::Int, Vpr::Matrix{Float64}, terms::NamedTuple,
                   vTol::Float64; verbose::Bool=false)

    beta   = terms.beta
    sigma  = terms.sigma
    phi    = terms.phi
    lgrid  = terms.lgrid
    agrid  = terms.agrid
    pil    = terms.pil
    captax = terms.captax
    lamval = terms.lamval
    tau    = terms.tau
    r      = terms.r
    w      = terms.w

    V  = zeros(nl, na)
    G  = zeros(nl, na)
    TV = zeros(nl, na)
    TG = zeros(nl, na)
    V0 = zeros(nl, na)

    # warm-start
    scale = 0.25
    for ia in 1:na
        kval = agrid[ia]
        for il in 1:nl
            yval = scale * (1.0 + r * (1.0 - captax[il])) * kval +
                   w * lgrid[il] - r * phi
            V[il, ia] = log(max(1e-10, yval))
        end
    end

    # expected values from continuation Vpr
    for ia in 1:na
        for il in 1:nl
            V0[il, ia] = dot(pil[il, :], Vpr[:, ia])
        end
    end

    # GSS over each (il, ia)
    for ia in 1:na
        for il in 1:nl
            y = w * lgrid[il] + (1.0 + r * agrid[ia])
            Cvals = V0[il, :]
            v, g, _ = gss(Cvals, y, beta, sigma, agrid, vTol * 1e-3)
            V[il, ia] = v
            G[il, ia] = g
        end
    end

    if verbose
        d    = supnorm(TV, V, 2)
        kdst = supnorm(TG, G, 2)
        @printf("\n\tBacksolve: ||TV - V|| = %4.6f\t||TG - G|| = %4.6f\n", d, kdst)
    end

    return V, G
end

end