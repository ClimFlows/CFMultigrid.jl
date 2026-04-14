using CFMultigrid
using Test

function test_solve(shape,smoother,kwargs,OPE,BC,verbose)
    # nx,ny,nz = 128,128,8
    #nx,ny,nz = 512,512,1
    #nx,ny,nz = 32,32,1
    nx,ny,nz = shape
    nhalo = 1

    #levels = setup_levels(nx,ny,nz)

    # smoother = Jacobi(0.85)
    # smoother = LinearRelaxation(nz,Float64)
    param = Param(1,2,20,10,1e-9,smoother)

    # kwargs = (;cxx=1,cyy=1,czz=100)
    # kwargs = (;)

    mg = setup(shape,nhalo,OPE,BC;kwargs...)

    (;b) = mg[1]

    b[1+div(nx,3),1+div(ny,3),1] = 1
    b[nx-div(nx,3),ny-div(ny,3),1] = -1

    solve!(mg,param,VCYCLE;verbose=verbose)
end


function test_default()
    res, ite = test_solve((128,128,8),Jacobi(0.85),(;),Poisson,NEUMANN, false)
    #res, ite = test_solve((128,128,8),LineRelaxation(8,Float64),(;cxx=1,cyy=1,czz=100),PoissonNonUniform,NEUMANN, true)
    #res, ite = test_solve((128,128,8),Jacobi(0.85),(;cxx=1,cyy=1,czz=100),PoissonNonUniform,NEUMANN, true)
     #res, ite = test_solve((128,128,16),Jacobi(0.85),(;),Poisson,DIRICHLET, true)

    @test isapprox(res, 1.01e-10, rtol = 0.05) & (ite==5)
end

function test_solve()
    exps = [
        [((64,64,64),LineRelaxation(64,Float64),(;),Poisson,NEUMANN),(3.25e-11, 5)],
        [((64,64,64),Jacobi(0.85),(;),Poisson,DIRICHLET),(1.11e-10, 6)],
        [((64,64,64),LineRelaxation(64,Float64),(;),Poisson,DIRICHLET), (1.47e-10, 5)],
        [((128,128,8),LineRelaxation(8,Float64),(;cxx=1,cyy=1,czz=100),PoissonNonUniform,NEUMANN),(1.45e-10,4)]
    ]
    for (a,out) in  exps
        res, ite = test_solve(a..., false)
        @test (isapprox(res,out[1], rtol=0.05) & (ite==out[2]))
    end
end

function test_helmholtz()
    exps = [
        [((128,128,8),LineRelaxation(8,Float64),(;cxx=1,cyy=1,czz=100,coef=1.0),Helmholtz,NEUMANN), (2.8e-10,5)]
    ]
    for (a,out) in  exps
        res, ite = test_solve(a..., true)
        @test (isapprox(res,out[1], rtol=0.05) & (ite==out[2]))
    end
end
