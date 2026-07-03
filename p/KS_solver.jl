
# making sure necessary packages are installed
using Pkg
Pkg.activate(".")   # so you're using the same packages as me
Pkg.instantiate()

using JLD2
using Random
using Printf
using LinearAlgebra: dot
using StatsBase: countmap
using Dates

include("src/ModelTypes.jl")
@printf("Loaded ModelTypes.jl\n")

include("src/Compute.jl")
@printf("Loaded Compute.jl\n")

include("src/ModelFunctions.jl")
@printf("Loaded ModelFunctions.jl\n")

include("src/EGM.jl")
@printf("Loaded EGM.jl\n")

include("src/DistrTools.jl")
@printf("Loaded DistrTools.jl\n")

include("src/Solvers.jl")
@printf("Loaded Solvers.jl\n")

include("src/SteadyState.jl")
@printf("Loaded Steady States.jl\n")

include("src/Predict.jl")
@printf("Loaded Predict.jl\n")

using .ModelTypes
using .ModelFunctions
using .Compute
using .EGM
using .DistrTools
using .Solvers
using .SteadyState
using .Predict

# model parameters
const α::Float64 = 0.36;
const β::Float64 = 0.96;
const δ::Float64 = 0.06;
const σ::Float64 = 2;
const ϕ::Float64 = 0;

# grid sizes and parameters
const na::Int64 = 100; 
const nl::Int64 = 7;
const nz::Int64 = 5;
const nk::Int64 = 25;

const a_l::Float64 = 0;
const a_h::Float64 = 100;

# calibrations for idiosyncratic income (vibes)
const μ_l::Float64 = 0; 
const ρ_l::Float64 = .9;
const σ_l::Float64 = .2;

# calibrations for aggregate TFP (Khan and Thomas 2013 would have 
# ρ_z = 0.909; I am setting it lower for two dimensional z with
# some variance for now.)
const μ_z::Float64 = 0;
const ρ_z::Float64 = .909;
const σ_z::Float64 = 0.014;

# kgrid
const kL::Float64 = 6; const kH::Float64 = 28; # this could be informed by steady states, but for now we do it this way
Kgrid = collect(range(kL, kH, length = nk));

grid_range = 2.575;

π_l, lgrid = getTauchen(nl,  μ_l, σ_l, ρ_l, grid_range);
π_z, zgrid= getTauchen(nz,  μ_z, σ_z, ρ_z, grid_range);

stationary_l = stationary(π_l);
const lagg::Float64 = dot(stationary_l, lgrid);
 
agrid = logspace(a_l, a_h, na);
amu = collect(range(a_l, a_h, length=na*10));

const np::Int64 = 10; # number of policies
const pol_l::Float64 = 0;
const pol_h::Float64 = .2;

τ_grid = range(pol_l, pol_h, length = np);
η_grid = range(pol_l, pol_h, length = np);
policy_grid = [(η, τ) for η in η_grid, τ in τ_grid];
captax = repeat([0.0], outer = nl);

const NT = 5000; #three thousand periods for sampling
const rnseed = 1234567;

const zt = simz(NT, nz, rnseed, π_z);

const params = ModelParams(α, β, δ, σ, ϕ, agrid, 
		lgrid, zgrid, π_l, π_z, amu, Kgrid);

const Kfore_start = [0.102898 0.946118; 
			0.112504 .944115;
			0.121485 0.942448;
			0.131148 0.940548;
			0.139579 0.939373]

const dTol = 1e-3;
const vTol = 1e-6;

results = Dict{Tuple{Float64,Float64}, NamedTuple}()
datestr = Dates.format(now(), "yyyy-mm-dd_HHMM")

for idx in eachindex(policy_grid)[1:4]
	η, τ = policy_grid[idx]

	@printf("Solving (eta, tau) = (%4.2f, %4.2f)--------------\n", η, τ)

	# set up local forecast to be updated within each thread
	# pulling from the forecast for the non-policy distortion version

	r_vals = zeros(nk, nz); w_vals = zeros(nk, nz); λ_vals = zeros(nk, nz); 

	for ik = 1:nk, iz = 1:nz
		r_vals[ik, iz] = calcr(α, δ, Kgrid[ik], η, zgrid[iz])
		w_vals[ik, iz] = calcw(α, Kgrid[ik], η, zgrid[iz])
		denom = dot((w_vals[ik, iz] .* lgrid).^(1 - τ), stationary(π_l))
		tot_inc = w_vals[ik, iz] * dot(lgrid, stationary(π_l))
		λ_vals[ik, iz] = tot_inc / denom
	end

	local policies = ProposedPolicies(η, τ, captax);

	local prices = ImpliedRegimeParams_KS(λ_vals, r_vals, w_vals)

	# init V and friends:
	V0 = zeros(nk,nz,nl,na); V  = zeros(nk,nz,nl,na);
	G0 = zeros(nk,nz,nl,na); G  = zeros(nk,nz,nl,na);  
	C  = zeros(nk,nz,nl,na);

	for ik = 1:nk, iz = 1:nz, il = 1:nl, ia = 1:na
		kval = agrid[ia];
		yval = (1 + r_vals[ik, iz]*(1-captax[il]))*kval + w_vals[ik, iz]*lgrid[il] - r_vals[ik, iz]*ϕ;
		ymin = max(1e-10, yval);
		V0[ik, iz, il, ia] = log(ymin);
		G0[ik, iz, il, ia] = agrid[ia];
	end
	
	@printf("K range for (%4.2f, %4.2f): %2.4f, %2.4f\n", 
		η, τ, minimum(Kt), maximum(Kt))
	Kfore_out, Kt = run_KS(V, V0, G, G0, C, params, policies, prices,
                zt, Kfore_start, vTol, dTol, verbose = true)
	
	# storing each stationary equilibrium
	results[(η, τ)] = (Kfore = Kfore_out, V = V, G = G, Kt = Kt,
		policies = policies)

end

@save "KS_Solves_$(datestr).jld2" results params zt
println("Done. Saved to KS_Solves_$(datestr).jld2")
