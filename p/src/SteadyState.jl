module SteadyState

export solveHousehold
import ..ModelTypes: ModelParams, ImpliedRegimeParams, ProposedPolicies
import ..Compute: supnorm, stationary
import ..EGM: solve
import ..ModelFunctions: calcr, calcw, tax
import ..Printf: @printf
import LinearAlgebra: dot

# solve for steady state distribution and value function, given a guess for the steady state capital stock. 
# Returns (V, G, C).
function solveHousehold(model::ModelParams, policies::ProposedPolicies, kval, vTol)

	π_l = model.pil;
	π_z = model.piz;

	lgrid = model.lgrid;
	agrid = model.agrid;
	zgrid = model.zgrid;

	# turns out I had this all wrong, so inputting the updated part
	# π_z is nz×nz over (yesterday's z, today's z). We track the *pair*
	# (z₋₁, z₀) as the aggregate state, so each household state carries a
	# transition index it ∈ 1:nt, nt = nz².
	#
	# CI = CartesianIndices(π_z):  CI[it] = (i, j) = (yesterday z_i, today z_j)
	#   → row    i = CI[it][1] = z₋₁  (yesterday)
	#   → column j = CI[it][2] = z₀   (today)   ← prices/flow payoff use THIS
	#
	# LI = LinearIndices(π_z): column-major, first index fastest.
	#   (1,1)→1  (2,1)→2  (1,2)→3  (2,2)→4   for nz=2
	#   i.e. LI[i,j] = i + (j-1)*nz
	#
	# Rolling forward: today's z_j becomes tomorrow's yesterday, so the
	# next pair is (j, k) = LI[j, k_next]  — indexed by the COLUMN, today.

	CI = CartesianIndices(π_z);
	LI = LinearIndices(π_z);
	nt = length(LI); # number of potential transitions

	# init value function
	nl = length(lgrid); na = length(agrid); nz = length(zgrid);
	V = zeros(nt, nl, na); V0 = zeros(nt, nl, na);
	EV = zeros(nt, nl, na);
	G = zeros(nt, nl, na); G0 = zeros(nt, nl, na);
	C = zeros(nt, nl, na);

	α = model.α; δ = model.δ; ϕ = model.ϕ;
	η = policies.η; τ = policies.τ; captax = policies.captax;
	# init price vectors
	r_val = calcr(α, δ, kval, η, zgrid); w_val = calcw(α, kval, η, zgrid);

	# we have to get the value of lambda such that taxation 
	# is redistributing everything. aka budget balance
	
	tot_inc = w_val.* dot(lgrid,stationary(π_l));
	denom = zeros(nz);
	for iz in 1:nz
		denom[iz] = dot((w_val[iz] * lgrid) .^ (1 - policies.τ),stationary(π_l));
	end
	λ = tot_inc./denom;

	prices = ImpliedRegimeParams(λ, r_val, w_val, captax);

	for it = 1:nt
		for il = 1:nl
			for ia = 1:na
				kval = agrid[ia];
				yval = (1 + r_val[1]*(1-captax[il]))*kval + w_val[1]*lgrid[il] - r_val[1]*ϕ;
				ymin = max(1e-10, yval);
				V0[it, il, ia] = log(ymin);
				G0[it, il, ia] = agrid[ia];
			end
		end
	end

	print("Solving Household Problem...\n")

	vdist::Float64 = 10;
	iter_ct = 1;
	while vdist > vTol

		for it in 1:nt
			today = CI[it][2] # <-- takes ONLY the column = today's z for now; need to update this once I do take previous z into account
			for il in 1:nl
				for ia in 1:na
					ev = 0.0
					for iz_next in 1:nz
						it_next = LI[today, iz_next]
						ev += π_z[today, iz_next] * dot(π_l[il, :], V0[it_next, :, ia])
					end
					EV[it, il, ia] = ev
				end
			end
		end

		V, G, C = solve(EV, model, policies, prices, CI)

		vdist = max(supnorm(V, V0, 3), supnorm(G, G0, 3))
		if iter_ct % 200 == 0
			@printf("\tIteration %i: ||V - V0|| = %1.6f, ||G - G0|| = %1.6f, dist = %1.6f\n", 
			iter_ct, supnorm(V,V0,3), supnorm(G,G0,3), vdist)
		end
		
		V0 .= V
		G0 .= G
	
		iter_ct += 1;
	end

	@printf("\tIteration %i: ||V - V0|| = %1.6f, ||G - G0|| = %1.6f, dist = %1.6f\n", 
				iter_ct, supnorm(V,V0,3), supnorm(G,G0,3), vdist)
	
	return V, G, C, CI, LI

end

end