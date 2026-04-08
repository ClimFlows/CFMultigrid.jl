using CFMultigrid
using Test

include("tests.jl")

@testset verbose=true "CFMultigrid" begin
    @testset "solve" begin
        test_default()
        test_solve()
    end
end
println("tests done")
