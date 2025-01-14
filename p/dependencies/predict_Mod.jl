module predict

include("HH_Mod.jl")
include("aiyagari_Mod.jl")
using .HH
using .aiyagari
using LinearAlgebra

export simz, transition, perfectForesight

function simz(N_T, Nz, rnseed, P_mat)
	zsim = zeros(1, N_T)
	izsim = zsim # indices

	Random.seed!(rnseed)
	efsim = rand(1, N_T)

	# start at median
	izsim[1] = div(Nz, 2)

	# build CDF
	cumP_mat = cumsum(P_mat, dims=2)

	# now, we simulate
	for t = 1:N_T-1
		csumvec = cumP_mat[izsim[t], 1:Nz]
		condmet = efsim[t+1] .<= csumvec
		izsim[t+1] = findfirst(condmet)
	end

	return izsim
end

function transition(r0, r1, lt, terms, lambda, dTol)
	verbose = false

	alpha = terms["alpha"]
	delta = terms["delta"]

	agrid = terms["agrid"]
	lgrid = terms["lgrid"]
	pil = terms["pil"]

	nl = length(lgrid)
	na = length(agrid)
	amu = range(minimum(agrid), maximum(agrid), length=na*10)

	T = length(lt)
	rt = range(only(r0), only(r1), T)
	rt = vcat(rt)
	wt = aiyagari.getW.(rt, alpha, delta)

	# initializing guess for K
	impliedK = rt .+ delta 
	impliedK = (impliedK)./ alpha
	impliedK = impliedK .^ (1 / (alpha - 1))
	impliedK = impliedK .* lt
	Kguess = impliedK

	# initializing storage arrays
	Garray = Vector{Any}(undef, T)
	Varray = Vector{Any}(undef, T)
	Farray = Vector{Any}(undef, T)

	# solving final period value function
	terms["r"] = rt[T]
	terms["w"] = wt[T]
	V, G, _ = HH.solve(nl, na, terms, dTol, verbose)
	Garray[T] = G
	Varray[T] = V
	adistr, Kagg = HH.getDist(G, amu, agrid, pil, true)
	Kguess[T] = Kagg
	Farray[T] = adistr

	# steady state in period 1
	terms["r"] = rt[1]
	terms["w"] = wt[1]
	V, G, _ = HH.solve(nl, na, terms, dTol, verbose)
	Garray[1] = G
	Varray[1] = V
	adistr, Kagg = HH.getDist(G, amu, agrid, pil, true)
	Kguess[1] = Kagg
	Farray[1] = adistr

	DIST = 10
	iter_ct = 1

	while DIST > dTol
		println("_Iteration $iter_ct")
		for t = T-1:-1:2
			if t % 50 == 0
				println("T = $t...")
			end

			Vpr = Varray[t+1]
			localterms = terms
			localterms["r"] = rt[t]
			localterms["w"] = wt[t]
			V, G, _ = HH.backsolve(nl, na, Vpr, localterms, dTol/10, verbose)

			Varray[t] = V
			Garray[t] = G
		end

		for t = 2:T-1
			if t % 50 == 0
				println("T = $t...")
			end
			Farray[t], Kguess[t] = HH.transitDistr(Garray[t], Farray[t-1], amu, agrid, pil)
		end

		DIST = LinearAlgebra.norm(Kguess .- impliedK, 1)
		println("||K - K'|| = $DIST")

		print(size(Kguess))
		print(size(lt))
		rguess = alpha .* (Kguess ./ lt) .^ (alpha - 1) - delta

		rt[2:T] = (1 - lambda) .* rt[2:T] + lambda .* rguess[2:T]
		wt = aiyagari.getW.(rt, alpha, delta)
		impliedK = (rt + delta) / alpha
		impliedK = impliedK .^ (1 / (alpha - 1))
		impliedK = impliedK .* lt
		iter_ct += 1
	end
end

function perfectForesight(zt, kst, k0, lambda, params, dTol)
	alpha = params[1]
	beta = params[2]
	delta = params[3]
	sigma = params[4]
	L = params[5]

	ist = delta * kst
	cst = kst^alpha - ist

	T = length(zt)
	zt = zt[1:T-1]
	kguess = fill(kst, T-1)
	kguess[1] = k0 * kst

	kprguess = fill(0, T-1)

	dist = maximum(abs.(kguess - kprguess))

	iter_ct = 1
	while dist > dTol
		yt = zt .* kguess .^ alpha .* (L^(-alpha))
		rt = alpha .* zt .* (kguess .^ (alpha - 1)) .* (L^(1 - alpha)) - delta

		ct = fill(cst, T-1)
		for i = T-2:-1:1
			if sigma == 1
				ct[i] = ct[i+1] * ((beta * (1 + rt[i+1]))^(-1))
			else
				ct[i] = ct[i+1] * ((beta * (1 + rt[i+1]))^(-1 / sigma))
			end
		end
		ct[T-1] = cst

		kt = fill(kst, T)
		kt[1] = k0 * kst
		for i = 1:T-1
			kt[i+1] = yt[i] - ct[i] + (1 - delta) * kt[i]
		end

		kprguess = lambda * kguess + (1 - lambda) * kt[1:T-1]

		dist = maximum(abs.(kprguess - kguess))
		println("Iteration $iter_ct: ||K' - K|| = $dist")

		iter_ct += 1
		kguess = kprguess
	end
end

end # module
