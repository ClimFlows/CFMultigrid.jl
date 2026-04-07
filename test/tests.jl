using CFMultigrid

function test_solve(shape,smoother,kwargs,OPE,BC)
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

    solve!(mg,param,VCYCLE)
end


function test_default()
    res, ite = test_solve((128,128,8),Jacobi(0.85),(;),Poisson,NEUMANN)
    return isapprox(res, 9.45e-11, atol = 1e-12) & (ite==5)
end
