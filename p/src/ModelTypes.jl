# ModelTypes.jl
module ModelTypes

export ModelParams, ImpliedRegimeParams, ProposedPolicies

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
    pil::Matrix{Float64}
    piz::Matrix{Float64}

    # distribution grid
    amu::Vector{Float64}
end

struct ImpliedRegimeParams
    # regime-specific parameters, which are fixed in each regime but can differ across regimes  
    λ::Vector{Float64}

    # resultant prices and wages, which depend on the regime
    r::Vector{Float64}
    w::Vector{Float64}

    # this just has to be dropped in again
    captax::Vector{Float64} 
end

struct ProposedPolicies
    # migration and progressivity
    η::Float64
    τ::Float64

    #capital income tax
    captax::Vector{Float64}
end

end