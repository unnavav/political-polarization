# imports everything in the right order, so that we can use the same names for structs and functions across files without worrying about circular dependencies

module PopulismModel

using Printf: @printf

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

# re-export so the notebook can access everything
using .ModelTypes
using .ModelFunctions
using .Compute	
using .EGM
using .DistrTools
using .Solvers

end