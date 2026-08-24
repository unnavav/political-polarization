# ModelTypes.jl
module ModelTypes

export ModelParams, ImpliedRegimeParams, ImpliedRegimeParams_KS, ProposedPolicies

struct ModelParams
    # household preference parameters
    α::Float64
    β::Float64
    δ::Float64
    σ::Float64

    # borrowing constraint
    ϕ::Float64

    # grids and transition probabilities
    agrid::Vector{Float64}
    lgrid::Vector{Float64}
    zgrid::Vector{Float64}
    Kgrid::Vector{Float64}
    Θgrid::Vector{Float64}

    amu::Vector{Float64}
    kernel::Vector{Float64}
    
    π_l::Matrix{Float64}
    π_z::Matrix{Float64}
    π_Θ::Matrix{Float64}
end

struct ImpliedRegimeParams_KS{N}
    # regime-specific parameters, which are fixed in each regime but can differ across regimes  
    λ::Array{Float64,N}

    # resultant prices and wages, which depend on the regime
    r::Array{Float64,N}
    w::Array{Float64,N}
end

struct ImpliedRegimeParams
    # regime-specific parameters, which are fixed in each regime but can differ across regimes  
    λ::Vector{Float64}

    # resultant prices and wages, which depend on the regime
    r::Vector{Float64}
    w::Vector{Float64}
end

struct ProposedPolicies
    # migration and progressivity
    η::Float64
    τ::Float64

    #capital income tax
    captax::Vector{Float64}
end

end