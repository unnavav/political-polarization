module compute

include("aiyagari_Mod.jl")
using .aiyagari

using StaticArrays

export gss, linterpolate, getKgrid, logspace, getTauchen, weight, cubicSpline, getSplineVal, dist, condense

function gss(Cvals, params, searchgrid, prec)
	ni = length(searchgrid)
	r = (3 - sqrt(5)) / 2

	y = params[1]

	a = searchgrid[1]
	b = min(searchgrid[ni], y) # Ensure c > 0
	c = (1 - r) * a + r * b
	d = r * a + (1 - r) * b

	vc = linterpolate(Cvals, searchgrid, c)
	fc = -aiyagari.util(c, params, vc)

	vd = linterpolate(Cvals, searchgrid, d)
	fd = -aiyagari.util(d, params, vd)
	iter_ct = 1

	v = [c fc; d fd]

	while abs(a - b) > prec
		if fc > fd
			a = c
			c = d
			d = r * a + (1 - r) * b

			fc = fd
			vd = linterpolate(Cvals, searchgrid, d)
			fd = -aiyagari.util(d, params, vd)
			v = [v; d fd]
		else
			b = d
			d = c
			c = (1 - r) * a + r * b

			fd = fc
			vc = linterpolate(Cvals, searchgrid, c)
			fc = -aiyagari.util(c, params, vc)
			v = [v; c fc]
		end

		iter_ct += 1
	end

	aval = c
	vc = linterpolate(Cvals, searchgrid, c)
	res = aiyagari.util(c, params, vc)

	return res, aval, v
end

function linterpolate(Vvec, searchgrid, vi)
	ni = length(Vvec)
	if vi <= searchgrid[1]
		v = Vvec[1]
	elseif vi >= searchgrid[ni]
		v = Vvec[ni]
	else
		a = collect(searchgrid .< vi)
		il = a |> findlast
		il = il[2]
		wl = (searchgrid[il + 1] - vi) / (searchgrid[il + 1] - searchgrid[il])
		v = wl * Vvec[il] + (1 - wl) * Vvec[il + 1]
	end

	return v
end

function getKgrid(nk, kl, kh)
	kgrid = logspace(log(kl - kl + 1) / log(10.0), log(kh - kl + 1) / log(10.0), nk)'
	kgrid += ones(size(kgrid)) * (kl - 1)
	kgrid = kgrid'
	
	return kgrid
end

function logspace(l, h, nx)
	grid = logspace(log(l + -1.0 * l + 1.0) / log(10.0), log(h + -1.0 * l + 1.0) / log(10.0), nx)'
	grid += ones(size(grid)) * (l - 1.0)

	grid = SVector{nx, Float64}(grid)
	return grid
end

function getTauchen(Nz, mu, sigma, rho)
	s = 2.575
	sigma_x = ((sigma^2) / (1 - rho^2))^(0.5)

	x_1 = mu - s * sigma_x
	x_Nz = mu + s * sigma_x

	x_grid = range(x_1, stop=x_Nz, length=Nz)
	z_grid = exp.(x_grid)

	w = x_grid[2] - x_grid[1]

	P_mat = MMatrix{Nz, Nz, Float64}zeros(Nz, Nz)

	for r in 1:Nz
		x_curr = x_grid[r] * rho

		P_mat[r, 1] = normcdf(x_grid[1] - x_curr + w / 2, mu, sigma)
		P_mat[r, Nz] = 1 - normcdf(x_grid[Nz] - x_curr - w / 2, mu, sigma)

		for c in 2:Nz-1
			upper = normcdf(x_grid[c] - x_curr + w / 2, mu, sigma)
			lower = normcdf(x_grid[c] - x_curr - w / 2, mu, sigma)
			P_mat[r, c] = upper - lower
		end
	end

	P_mat = SMatrix{Nz, Nz, Float64}(P_mat)
	z_grid = SVector{Nz, Float64}(z_grid)

	return P_mat, z_grid
end

function weight(pim, f)
	if f >= maximum(pim)
		ix = length(pim) - 1
		we = 0
	elseif f < minimum(pim)
		ix = 1
		we = 1
	else
		ix = sum(pim .<= f)
		we = (pim[ix + 1] - f) / (pim[ix + 1] - pim[ix])
	end

	return ix, we
end

function cubicSpline(l, h, r, VK)
	lb = log(l) / log(10)
	ub = log(h) / log(10)

	val_grid = logspace(lb, ub, r + 2)

	ffunc(ki) = VK[1, ki]

	ti = val_grid
	fti = zeros(r + 2)
	dti = zeros(r + 2)
	ttf = zeros(r + 2)

	fti[1] = ffunc(1)

	for i in 2:r + 2
		t = val_grid[i]
		t_1 = val_grid[i - 1]

		dti[i] = t - t_1
		fti[i] = ffunc(i)
		ttf[i] = (fti[i] - fti[i - 1]) / dti[i]
	end

	upper_diag = zeros(r + 1)
	lower_diag = zeros(r + 1)
	principal = zeros(r + 1)

	for i in 2:r
		upper_diag[i] = dti[i]
		lower_diag[i] = dti[i + 2]
	end

	upper_diag[2] += (dti[2]^2) / dti[3]
	lower_diag[r] += (dti[r]^2) / dti[r - 1]

	omega_1 = dti[3] - (dti[2]^2) / dti[3]
	omega_r = dti[r + 1] - (dti[r + 2]^2) / dti[r + 1]

	principal[1] = 2 * (dti[2] + dti[3]) - omega_1
	principal[r + 1] = 2 * (dti[r + 1] + dti[r + 2]) - omega_r

	for i in 2:r
		principal[i] = 2 * (dti[i + 1] + dti[i + 2])
	end

	f = zeros(r)

	f[1] = 3 * (dti[3] * ttf[2] + dti[2] * ttf[3]) -
		2 * (dti[3] * ttf[2] - (dti[2]^2) / dti[3] * ttf[3])

	f[r] = 3 * (dti[r + 2] * ttf[r + 1] + dti[r + 1] * ttf[r + 2]) -
		2 * (dti[r + 1] * ttf[r + 2] - (dti[r + 2]^2) / dti[r + 1] * ttf[r + 1])

	for i in 2:r
		f[i] = 3 * (dti[i + 2] * ttf[i + 1] + dti[i + 1] * ttf[i + 2])
	end

	T = spdiagm(-1 => circshift(lower_diag, -1), 0 => principal, 1 => circshift(upper_diag, 1))
	T[r, r] = 2 * (dti[r + 1] + dti[r + 2]) - omega_r

	s = T \ f

	s_0 = 2 * ttf[2] - ((dti[2] / dti[3])^2) * ttf[3] -
		(1 - (dti[2] / dti[3]^2)) * s[1] +
		(dti[2] / dti[3]^2) * s[2]

	s_r1 = 2 * (ttf[r + 2] - ((dti[r + 2] / dti[r + 1]))^2 * ttf[r + 1]) -
		(1 - (dti[r + 2] / dti[r + 1])^2) * s[r] +
		(dti[r + 2] / dti[r + 1])^2 * s[r - 1]

	s_fin = [s_0; s; s_r1]

	C = zeros(r + 1, 4)

	for i in 1:r + 1
		C[i, 1] = fti[i]
		C[i, 2] = s_fin[i]
		C[i, 3] = 3 * ttf[i + 1] / dti[i + 1] -
			2 * s_fin[i] / dti[i + 1] - s_fin[i + 1] / dti[i + 1]
		C[i, 4] = (-2 * ttf[i + 1] +
			s_fin[i] + s_fin[i + 1]) / (dti[i + 1]^2)
	end

	return C
end

function getSplineVal(coeffs, lval, lgrid)
	il = sum(lgrid .<= lval)
	x = lgrid[il]
	if il == length(lgrid)
		il -= 1
	end

	d = lval - x

	coeffs = squeeze(coeffs)

	val = coeffs[1, il] +
		coeffs[2, il] * d +
		coeffs[3, il] * d^2 +
		coeffs[4, il] * d^3

	return val
end

function dist(M, N, nd)
	x = abs.(M .- N)
	for i in 1:nd
		x = max(x)
	end
	return x
end

function condense(adistr, amu, agrid)
	nl, nmu, np = size(adistr)
	na = length(agrid)
	distr = zeros(nl, na, np)
	for im in 1:nmu
		ix, we = weight(agrid, amu[im])

		for ip in 1:np
			distr[:, ix, ip] += we * adistr[:, im, ip]
			distr[:, ix + 1, ip] += (1 - we) * adistr[:, im, ip]
		end
	end
	return distr
end

end # module
