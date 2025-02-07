module HH

export solve, backsolve, getDist, transitDistr, map, util

include("egm_Mod.jl")
include("gov_Mod.jl")
include("compute_Mod.jl")
using .egm
using .gov
using .compute

using LinearAlgebra, Printf, StaticArrays, FLoops

function solve(nl, na, terms, vTol, verbose)
	beta = terms["beta"]
	sigma = terms["sigma"]
	phi = terms["phi"]
	lgrid = terms["lgrid"]
	agrid = terms["agrid"]
	pil = terms["pil"]
	g = terms["G"]
	captax = terms["captax"]
	lamval = terms["lamval"]
	tau = terms["tau"]

	r = terms["r"]
	w = terms["w"]

	V = zeros(nl, na)
	G = zeros(nl, na)
	TG = zeros(nl, na)
	TV = zeros(nl, na)
	V0 = zeros(nl, na)
	endoK = zeros(nl, na)
	C = zeros(nl, na)

	# set up V so that it doesn't start empty
	scale = 0.25
	@floop ThreadedEx() for ia = 1:na, il = 1:nl
		kval = agrid[ia]
		yval = scale * (1 + r * (1 - captax[il])) * kval + w * lgrid[il] - r * phi
		ymin = max(1e-10, yval)
		V[il, ia] = log(ymin)
	end

	# init expected vals
	@floop ThreadedEx() for ia = 1:na, il = 1:nl
		V0[il, ia] = dot(pil[il, :], V[:, ia])
	end

	dist = 1e5
	kdist = dist
	iter_ct = 1

	# if verbose
		println("Begin EGM.")
	# end

	while dist > vTol

		# endogrid choice
		@floop ThreadedEx() for ia = 1:na, il = 1:nl
			kpr = agrid[ia]
			l = lgrid[il]
			D = egm.solveD(V0[il, :], ia, agrid)

			endoK[il, ia] = ((beta * D)^(-1) + kpr - w * l +
				gov.tax(w * l, lamval, tau) - g[1] +
				r * (1 - captax[il]) * phi) / (1 + r * (1 - captax[il]))
		end

		# linterpolate decision rule
		@floop ThreadedEx() for ia = 1:na, il = 1:nl
			k = agrid[ia]
			lkvals = endoK[il, :]

			if k < lkvals[1]
				G[il, ia] = agrid[1]
			else
				ix, we = compute.weight(lkvals, k)
				kpr = we * agrid[ix] + (1 - we) * agrid[ix + 1]
				G[il, ia] = max(agrid[1], kpr)
			end
		end

		@floop ThreadedEx() for ia = 1:na, il = 1:nl
			l = lgrid[il]

			tax_labor = gov.tax(w * l, lamval, tau)
			c = (1 + r * (1 - captax[il])) * agrid[ia] + w * l - tax_labor +
				g - G[il, ia] - r * (1 - captax[il]) * phi
			C[il, ia] = max(1e-10, c)

			ix, we = compute.weight(agrid, G[il, ia])
			ev = we * V0[il, ix] + (1 - we) * V0[il, ix + 1]
			TV[il, ia] = log(c) + beta * ev
		end

		dist = maximum(maximum(abs.(TV .- V)))
		kdist = norm(TG .- G, 2)

		iter_ct += 1

		V .= TV
		TG .= G

		@floop ThreadedEx() for ia = 1:na, il = 1:nl
			V0[il, ia] = dot(pil[il, :], V[:, ia])
		end

	end

	# if verbose
		println(@sprintf("\n\tIteration %i: ||TV - V|| = %1.8f\t||TG - G|| = %1.8f", iter_ct, dist, kdist))
	# end
	return V, G, V0

end

function backsolve(nl, na, Vpr, terms, vTol, verbose)
	beta = terms["beta"]
	sigma = terms["sigma"]
	phi = terms["phi"]
	lgrid = terms["lgrid"]
	agrid = terms["agrid"]
	pil = terms["pil"]
	g = terms["G"]
	captax = terms["captax"]
	lamval = terms["lamval"]
	tau = terms["tau"]

	r = terms["r"]
	w = terms["w"]

	V = zeros(nl, na)
	G = zeros(nl, na)
	TG = zeros(nl, na)
	TV = zeros(nl, na)
	V0 = zeros(nl, na)

	# set up V so that it doesn't start empty
	scale = 0.25
	@floop ThreadedEx() for ia = 1:na, il = 1:nl
		kval = agrid[ia]
		yval = scale * (1 + r * (1 - captax[il])) * kval + w * lgrid[il] - r * phi
		ymin = max(1e-10, yval)
		V[il, ia] = log(ymin)
	end

	# expected vals
	@floop ThreadedEx() for ia = 1:na, il = 1:nl
		V0[il, ia] = pil[il, :] ⋅ Vpr[:, ia]
	end

	# get best capital choice from V, EV using GSS
	@floop ThreadedEx() for ia = 1:na, il = 1:nl
		y = w * lgrid[il] + (1 + r * agrid[ia])
		Cvals = V0[il, :]
		params = [y, beta, sigma]

		v, g, ~ = compute.gss(Cvals, params, agrid, vTol / 1e3)

		V[il, ia] = v
		G[il, ia] = g
	end

	dist = maximum(maximum(abs.(TV .- V)))
	kdist = norm(TG .- G, 2)

	if verbose
		println("\n\tIteration $iter_ct: ||TV - V|| = $dist\t||TG - G|| = $kdist")
	end

	return V,G,V0
end

function getDist(G, amu, agrid, pil, verbose)
	nl, np = size(G)
	nmu = size(amu)[1]
	mu = zeros(nl, nmu)

	# nomenclature: ixagrid is indices for HH a for both parties.
	ixgrid = zeros(nl, nmu)
	wegrid = zeros(size(ixgrid))

	# We linearly interpolate the policy function on agrid to compute
	# afval on the distribution support.
	for im = 1:nmu
		kval = amu[im]
		ixk, wek = compute.weight(agrid, kval)

		for il = 1:nl

			kdval = G[il, ixk] * wek + G[il, ixk + 1] * (1.0 - wek)

			ix, we = compute.weight(amu, kdval)
			ixgrid[il, im] = ix
			wegrid[il, im] = we
		end
	end
	

	ixgrid = round.(Int,ixgrid)

	distance = 20
	iter_ct = 1

	# nomenclature: muA is the distribution over assets when A is
	# in power
	mu = ones(size(mu)) * (1 / (nmu * nl))

	while (distance > 1e-6 && iter_ct < 3000)

		mu1 = zeros(size(mu))

		for im = 1:nmu
			for il = 1:nl

				# need weighting for a and b HH under party A
				ix = ixgrid[il, im]
				we = wegrid[il, im]
				muval = mu[il, im]

				if muval > 0
					for jl = 1:nl

						# a households' movements
						if (ix < nmu)
							mu_val1 = pil[il, jl] * muval * we
							mu_val2 = pil[il, jl] * muval * (1.0 - we)
						else
							mu_val1 = pil[il, jl] * muval * we
							mu_val2 = 0
						end

						mu1[jl, ix] += mu_val1
						mu1[jl, ix + 1] += mu_val2
					end
				end
			end
		end

		distance = norm(mu1.-mu, 2)

		iter_ct += 1

		mu = mu1
	end

	s = "\n\t\tIteration $iter_ct: ||Tm - m|| = $distance\tsum = $(sum(sum(mu)))"
	if verbose
		println(s)
	end

	distrA2500 = dropdims(sum(mu, dims = 1), dims = 1)
	kagg = dot(collect(amu), distrA2500')
	# println(@sprintf("Kagg = %1.4f", kagg))
	# println(@sprintf("sum(mu) = %1.4f", sum(sum(mu))))

	return mu, kagg
end

function transitDistr(G, mu, amu, agrid, pil)
	nl, np = size(G)
	nmu = size(amu)[1]
	mu = ones(size(mu)) * (1 / (nmu * nl))

	# nomenclature: ixagrid is indices for HH a for both parties.
	ixgrid = zeros(nl, nmu)
	wegrid = zeros(size(ixgrid))

	# We linearly interpolate the policy function on agrid to compute
	# afval on the distribution support.
	for im = 1:nmu
		kval = amu[im]
		ixk, wek = compute.weight(agrid, kval)

		for il = 1:nl

			kdval = G[il, ixk] * wek + G[il, ixk + 1] * (1.0 - wek)

			ix, we = compute.weight(amu, kdval)
			ixgrid[il, im] = ix
			wegrid[il, im] = we
		end
	end

	ixgrid = round.(Int,ixgrid)
	mu1 = zeros(size(mu))

	for im = 1:nmu
		for il = 1:nl

			ix = ixgrid[il, im]
			we = wegrid[il, im]
			muval = mu[il, im]

			if muval > 0
				for jl = 1:nl

					# a households' movements
					if (ix < nmu)
						mu_val1 = pil[il, jl] * muval * we
						mu_val2 = pil[il, jl] * muval * (1.0 - we)
					else
						mu_val1 = pil[il, jl] * muval * we
						mu_val2 = 0
					end

					mu1[jl, ix] += mu_val1
					mu1[jl, ix + 1] += mu_val2
				end
			end
		end
	end

	
	distrA2500 = sum(mu, dims = 1)
	kagg = dot(collect(amu), distrA2500')
	# println(@sprintf("Kagg = %1.4f", kagg))
	# println(@sprintf("sum(mu) = %1.4f", sum(sum(mu))))

	return mu, kagg
end

function mapVotes(VOTES, amu, agrid, adistr, pctDem)
	nl, np, nm = size(adistr)
	ixgrid = zeros(nm)
	wegrid = ixgrid
	vdistr = zeros(nl, nm)

	# first get mapping
	for im = 1:nm
		a = amu[im]
		ix, we = compute.weight(agrid, a)

		ixgrid[im] = ix
		wegrid[im] = we
	end

	for im = 1:nm
		for ip = 1:np
			for il = 1:nl

				muval = adistr[il, ip, im]
				if muval > 0
					ix = ixgrid[im]
					we = wegrid[im]
					vdistr[il, im] = pctDem * VOTES[il, ix, 1] * we * muval +
						(1 - pctDem) * VOTES[il, ix, 2] * we * muval +
						pctDem * VOTES[il, ix + 1, 1] * (1 - we) * muval +
						(1 - pctDem) * VOTES[il, ix + 1, 2] * (1 - we) * muval
				end
			end
		end
	end

	majority = sum(sum(vdistr))

	if majority > 0.5
		winner = 1
	else
		winner = 0
	end
end

end  # module
