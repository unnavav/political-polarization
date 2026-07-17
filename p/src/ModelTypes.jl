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
    π_l::Matrix{Float64}
    π_z::Matrix{Float64}

    # distribution grid
    amu::Vector{Float64}
    # kgrid
    Kgrid::Vector{Float64}
end

struct ImpliedRegimeParams_KS
    # regime-specific parameters, which are fixed in each regime but can differ across regimes  
    λ::Matrix{Float64}

    # resultant prices and wages, which depend on the regime
    r::Matrix{Float64}
    w::Matrix{Float64}
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