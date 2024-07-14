using CFMultigrid
using Test

include("tests.jl")
include("test_RP.jl")
include("test_3d.jl")

@testset verbose=true "CFMultigrid" begin
    @testset "rectangular 2D" begin
        @test test_default()
        @test test_2d_centers()
        @test test_2d_vertices()
    end

    @testset "mask 2D" begin
        @test test_2d_msk_centers()
        @test test_2d_msk_vertices()
        @test test_2d_msk_vertices_convergence_vs_resolution()
        @test test_2d_msk_centers_convergence_vs_resolution()
        @test test_2d_msk_vertices_convergence_vs_rhs()
    end

    @testset "mask 3D" begin
        @test test_3d_centers()
        @test test_3d_vertices()

    end
end
println("tests done")
