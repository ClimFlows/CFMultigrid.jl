module CFMultigrid

export solve!, Param, setup_gmg, VCYCLE
export Jacobi, Gauss_Seidel, LineRelaxation
export Poisson, PoissonNonUniform
export NEUMANN, DIRICHLET
export restriction, prolongation

include("solver.jl")

end
