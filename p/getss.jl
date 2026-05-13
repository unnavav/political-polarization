# run if you don't have this package (I didn't), otherwise comment out
# import Pkg; Pkg.add("Printf"); Pkg.add("JLD2")
using JLD2

using Printf
using LinearAlgebra: dot
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

using .ModelTypes
using .ModelFunctions
using .Compute
using .EGM
using .DistrTools
using .Solvers
using .SteadyState

# model parameters
const α::Float64 = 0.36;
const β::Float64 = 0.96;
const δ::Float64 = 0.06;
const σ::Float64 = 2;
const ϕ::Float64 = 0;

# grid sizes and parameters
const na::Int64 = 100; 
const nl::Int64 = 7;
const nz::Int64 = 2;

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
const ρ_z::Float64 = .75;
const σ_z::Float64 = 0.014;


# need to back out σ^2_e given σ^2_l
σ_le = σ_l*sqrt(1 - ρ_l^2);
grid_range = 2.575;

π_l, lgrid = getTauchen(nl,  μ_l, σ_le, ρ_l, grid_range);

stationary_l = stationary(π_l);
const lagg::Float64 = dot(stationary_l, lgrid);

σ_ze = σ_z*sqrt(1 - ρ_z^2);
π_z, zgrid= getTauchen(nz,  μ_z, σ_ze, ρ_z, grid_range);

agrid = logspace(a_l, a_h, na);
amu = collect(range(a_l, a_h, length=na*10));

const params = ModelParams(α, β, δ, σ, ϕ, agrid, 
    lgrid, zgrid, π_l, π_z, amu);

const np::Int64 = 10; # number of policies
const pol_l::Float64 = 0;
const pol_h::Float64 = 1;

τ_grid = range(pol_l, pol_h, length = np)
η_grid = range(pol_l, pol_h, length = np)
policy_grid = [(η, τ) for η in η_grid for τ in τ_grid]

captax = repeat([0], outer = 7);

results = Dict{Tuple{Float64,Float64}, NamedTuple}()

datestr = Dates.format(now(), "yyyy-mm-dd_HHMM")
results_lock = ReentrantLock()

vTol::Float64 = 1e-6
kTol::Float64 = 1e-4

Threads.@threads for (η, τ) in policy_grid
	
	# all local to this thread
    local kl = 0.0
    local kh = 20.0
    local kval = 0.5 * (kl + kh)
    local kdist = 1e5
    local V = nothing
    local G = nothing
    local C = nothing
    local CI = nothing
    local LI = nothing
    local μ = nothing
    local K = 0.0

	policies = ProposedPolicies(η, τ, captax)
	
	@printf("Solving (eta, tau) = (%4.2f, %4.2f)--------------\n", η, τ)
	
	while kdist > kTol

		@printf("\n\nK guess: %4.8f\n", kval)
		
		V, G, C, CI, LI = solveHousehold(params, policies, kval, vTol);

		μ, K = getDistr(G, amu, agrid, π_l, π_z, CI, LI, ϕ, verbose = true, vTol = 1e-8);

		diff = K - kval
		adj = abs(diff) < .1 ? 0.5 : 0.7 # fancy if-then in one line

		if diff > 0
			#@printf("\n||K - kval|| = %4.5f. \tCapital too low.\n", abs(diff))
			kl = (1 - adj) * kval + adj * kl
		else
			#@printf("\n||K - kval|| = %4.5f. \tCapital too high.\n", abs(diff))
			kh = (1 - adj) * kval + adj * kh
		end

		kdist = abs(diff)
		kval = 0.5 * (kl + kh)

	end
	
	# storing each steady state
	lock(results_lock) do
        results[(η, τ)] = (V=V, G=G, C=C, μ=μ, K=kval, policies=policies)
    end

end

@save "steady_states_$(datestr).jld2" results params policy_grid
@printf("Done. Saved to steady_states_$(datestr).jld2\n")