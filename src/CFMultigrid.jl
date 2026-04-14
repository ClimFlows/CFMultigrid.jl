module CFMultigrid

include("axes.jl")
include("boundary_conditions.jl")
include("lowlevel_operators.jl")
include("grids.jl")
include("RPs.jl")
include("smoothers.jl")
include("poisson.jl")
include("poisson_nonuniform.jl")
include("helmholtz.jl")
include("highlevel_operators.jl")
include("setup.jl")
include("solver.jl")

export solve!, Param, setup, VCYCLE
export Jacobi, Gauss_Seidel, LineRelaxation
export Poisson, PoissonNonUniform, Helmholtz
export NEUMANN, DIRICHLET
export restriction, prolongation

end
