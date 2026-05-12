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

	# The cartesian index is the location, whereas the linear index is the one-dimensional 
	# version that determines whether we transition and from which state we do. So we have:
	#  - CI[i,j] = (i,j) is the cartesian index of the state (z_i, l_j)
	#  - LI[i,j] = k is the linear index of the state (z_i, l_j), 
	# 		which determines the transition probabilities. 
	# Because I'm interested in the combination of yesterday's TFP and today's TFP as a 
	# state variable, I need to use the cartesian index to determine the transition probabilities. 
	# The mapping for this is below. I am writing a long explanation here for my own sake. 
	#
	# This means (1,1) maps to 1 in the linear index, (1,2) maps to 2, (2,1) maps to 3, and 
	# (2,2) maps to 4, if nz = 2. Mapping back gives us:
	# 

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

	print("Solving Household Problem...")

	vdist::Float64 = 10;
	iter_ct = 1;
	while vdist > vTol

		for it in 1:nt
			today = CI[it][2]
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

		vdist = max(supnorm(V,V0,3), supnorm(G,G0,3));

		V0 .= .5*V + .5*V0; G0 .= .5*G + .5*G0;
		
		if iter_ct % 50 == 0
			@printf("Iteration %i: ||V - V0|| = %1.6f, ||G - G0|| = %1.6f, dist = %1.6f\n", 
			iter_ct, supnorm(V,V0,3), supnorm(G,G0,3), vdist)
		end
		iter_ct += 1;
	end

	@printf("Iteration %i: ||V - V0|| = %1.6f, ||G - G0|| = %1.6f, dist = %1.6f\n", 
				iter_ct, supnorm(V,V0,3), supnorm(G,G0,3), vdist)
	
	return V, G, C, CI, LI

end

end