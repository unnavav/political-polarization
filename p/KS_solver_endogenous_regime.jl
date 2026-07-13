# making sure necessary packages are installed

println("Adding Packages...\n")
#using Pkg
#Pkg.activate(".")   # so you're using the same packages as me
#rm("Manifest.toml", force=true)   # drop the 1.12-generated one
#Pkg.resolve()                      # rebuild for THIS Julia (1.11.7)
#Pkg.instantiate(; verbose=true)

using JLD2
using Random
using Printf
using LinearAlgebra
using StatsBase: countmap
using Dates

for file in [
    "src/ModelTypes.jl",
    "src/Compute.jl",
    "src/DistrTools.jl",
    "src/ModelFunctions.jl",
    "src/EGM.jl",
    "src/Solvers.jl",
    "src/SteadyState.jl",
    "src/Predict.jl"]
    include(file)
    @printf("Loaded %s\n", basename(file))
end

using .ModelTypes
using .Compute
using .EGM
using .DistrTools
using .Solvers
using .SteadyState
using .Predict
using .ModelFunctions

# ─── Frequency-invariant parameters ───
const α::Float64 = 0.36
const σ::Float64 = 2
const ϕ::Float64 = 0
const μ_l::Float64 = 0
const μ_z::Float64 = 0

# grid sizes and parameters
const na::Int64 = 100; 
const nl::Int64 = 15;
const nz::Int64 = 2;
const nk::Int64 = 37;

const a_l::Float64 = 0;
const a_h::Float64 = 100;

# policy grids
const np::Int64 = 10; # number of policies
const pol_l::Float64 = 0;
const pol_h::Float64 = .18;

τ_grid = range(pol_l, pol_h, length = np);
η_grid = range(pol_l, pol_h, length = np);
policy_grid = [(η, τ) for η in η_grid, τ in τ_grid];
captax = repeat([0.0], outer = nl);

# ─── Frequency switch: read from environment, default to "annual" ───
const freq = get(ENV, "FREQ", "quarterly")
display("Frequency set to: $freq")

if freq == "quarterly"
    const β::Float64   = 0.99
    const δ::Float64   = 0.025
    const ρ_l::Float64 = 0.9878      # STY persistence quarterly
    const σ_l::Float64 = 0.087       # STY innovation std, quarterly
    const ρ_z::Float64 = 0.976
    const σ_z::Float64 = 0.007
    
			
elseif freq == "annual"
    const β::Float64   = 0.96
    const δ::Float64   = 0.06
    const ρ_l::Float64 = 0.952     # Storesletten-Telmer-Yaron
    const σ_l::Float64 = 0.17       # = sqrt(0.061), STY persistent innovation std σ_η
    const ρ_z::Float64 = 0.909      # Khan-Thomas 2013
    const σ_z::Float64 = 0.014
    
else
    error("FREQ must be \"annual\" or \"quarterly\", got \"$freq\"")
end

const freq_label = freq   # use in output filename

# building fixed grids ───

grid_range_l = 3.5     # wide, for near-unit-root labor
grid_range_z = 2.575   # leave z as is

π_l, lgrid = getTauchen(nl,  μ_l, σ_l, ρ_l, grid_range_l);
# π_z, zgrid= getTauchen(nz,  μ_z, σ_z, ρ_z, grid_range_z);

# stationary_l = stationary(π_l);
# const lagg::Float64 = dot(stationary_l, lgrid);

zgrid = [0.98, 1.02]          # bad, good
π_z = [0.80  0.20;                 # from bad: 20% chance exit to good
       0.03125  0.96875]          # from good: ~3% chance enter recession
       
stationary_l = stationary(π_l);
const lagg::Float64 = dot(stationary_l, lgrid);
 
# more conservative estimates 
# const kH::Float64 = ((1.0/β - 1.0 + δ) / α)^(1.0/(α-1.0)) * (1.0 + maximum(η_grid)) * 1.5
# const kL::Float64 = max(1.0, ((1.0/β - 1.0 + δ) / α)^(1.0/(α-1.0)) * (1.0 + minimum(η_grid)) * 0.5)

const kH::Float64 = 50; const kL::Float64 = 42;
Kgrid = collect(range(kL, kH, length = nk));

Kfore_start = [0.081247 0.978478;
        0.090991 0.976496];

agrid = logspace(a_l, a_h, na);
amu = collect(range(a_l, a_h, length=na*10));

const NT = 5000; #three thousand periods for sampling
const rnseed = 1234567;

const zt = simz(NT, nz, rnseed, π_z);

const params = ModelParams(α, β, δ, σ, ϕ, agrid, 
		lgrid, zgrid, π_l, π_z, amu, Kgrid);

const dTol = 5e-3;
const vTol = 1e-5;

results = Dict{Tuple{Float64,Float64}, NamedTuple}()
datestr = Dates.format(now(), "yyyy-mm-dd_HHMM")

for idx in eachindex(policy_grid)
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
	
	Kfore_out, Kt = run_KS(V, V0, G, G0, C, params, policies, prices,
                zt, Kfore_start, vTol, dTol, λ_damp=0.5, verbose = true)
	
    flush(stdout)
    
    result_p = (Kfore = Kfore_out, V = V, G = G, Kt = Kt,
		 policies = policies);
    fname = @sprintf("policy_%.4f_%.4f.jld2", η, τ)
    @save "../d/ks/policy_results/$(fname)_$(freq_label)_nl$(nl)" result_p

    # storing each stationary equilibrium in mama file
	results[(η, τ)] = result_p

end

@save "../d/ks/KS_Solves_$(datestr)_nl$(nl).jld2" results params zt
println("Done. Saved to KS_Solves_$(freq_label)_$(datestr).jld2")