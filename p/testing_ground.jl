

using JLD2
using Random
using Printf
using LinearAlgebra: dot
using StatsBase: countmap
using Dates

for file in [
    "src/Compute.jl"]
    include(file)
    @printf("Loaded %s\n", basename(file))
end

using .Compute


# ─── Frequency-invariant parameters ───
const α::Float64 = 0.36
const σ::Float64 = 2
const ϕ::Float64 = 0
const μ_l::Float64 = 0
const μ_z::Float64 = 0

const β::Float64   = 0.99
const δ::Float64   = 0.025
const ρ_l::Float64 = 0.9878      # STY persistence quarterly
const σ_l::Float64 = 0.087       # STY innovation std, quarterly
const ρ_z::Float64 = 0.976
const σ_z::Float64 = 0.007

for gr in (1.0, 1.25, 1.5, 1.75, 2.0)
    πl, lg = getTauchen(7, μ_l, σ_l, ρ_l, gr)
    s = stationary(πl); mv = dot(s,(lg.-dot(s,lg)).^2); tv = σ_l^2/(1-ρ_l^2)
    @printf("gr=%.3f  ratio=%.3f\n", gr, mv/tv)
end