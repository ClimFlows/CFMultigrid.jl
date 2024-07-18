# using AnotherGeometricMultigrid
using CFMultigrid

function test_solve(;kwargs=Dict(),verbose=true)
    param = Param(;kwargs...)
    param.verbose = verbose
    mg = get_gmg(param)
    if verbose
        println(param)
        println(mg.levels)
    end
    grid = mg.levels[1].grid
    b = mg.data[1].b
    nh = grid.nh
    nx = grid.nx
    ny = grid.ny
    nz = grid.nz

    i0 = nh+6
    i1 = nh+4
    j0 = nh+div(ny,2)
    if mg.param.case == :threed
        k0 = nh+4
        b[i0,i0,k0] = 1.0
        b[i1,i1,k0] = -1.0
    else
        b[i0,j0,1] = 1.0
        b[i1,j0,1] = -1.0
    end

    res, nite = solve!(mg)
end


function test_default()
    res, ite = test_solve(;verbose=false)
    return isapprox(res, 9.52e-11, atol = 1e-13) & (ite==8)
end

function test_2d_centers()
    p = Dict([(:nx,512),
              (:ny,512),
              (:nz,1),
              (:case,:twod)])
    res, ite = test_solve(;verbose=false,kwargs=p)
    return isapprox(res, 3.210e-11, atol = 1e-13) & (ite==11)
end

function test_2d_vertices()
    p = Dict([(:nx,512),
              (:ny,512),
              (:nz,1),
              (:location, :vertices),
              (:case,:twod)])
    res, ite = test_solve(;verbose=false,kwargs=p)
    return isapprox(res, 1.976e-11, atol = 1e-13) & (ite==10)
end
