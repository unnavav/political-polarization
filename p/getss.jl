# GET SS . JL
# May 2026
# vaasavi, doing her best to translate to Julia 🫠

include("PopulismModel.jl")
using .PopulismModel

na = 100; nl = 7; nz = 5;

σ_l = .2;
ρ_l = .9;
# need to back out σ^2_e given σ^2_l
sigx = σ_l*sqrt(1 - ρ_l^2);
range = 2.575;
[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);


params = ModelParams(α=0.36, β=0.96, δ=0.08, σ=2.0, ϕ=0.0,
					 agrid=range(0.0, 50.0, length=na),
					 lgrid=lgrid,
					 zgrid=,
					 pil=, # placeholder
					 piz=)   # placeholder

Kmin = 3; Kmax = 20;



regime = RegimeParams(λ=0.5, τ=[0.2, 0.3], η=[0.1, 0.2], captax=repeat(0, nl),
				 r=, # placeholder
				 w=) # placeholder