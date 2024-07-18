include("test_3d.jl")

function test_3d_centers_rectangular()
    res, ite = run_rectangular(32,:centers)
    return isapprox(res, 2.6e-13, atol = 1e-14) & (ite==6)
end

function test_3d_vertices_rectangular()
    res, ite = run_rectangular(32,:vertices)
    return isapprox(res, 6.9e-13, atol = 1e-14) & (ite==6)
end


function run_rectangular(n,location;domain=:rectangular,config=:dipole)
    config = :dipole
    domain = :rectangular

    n = 32
    param = Param()
    param.nx = 8*n
    param.ny = n
    param.nz = n
    param.case=:threed
    param.location=location
    param.cycle=:v
    param.tol=1e-12
    param.ndeepest = 60
    param.npre = 2
    param.npost = 2
    param.omega = 0.88
    param.omega = 0.9
    param.maxite = 15
    param.verbose = true

    mg = get_gmg(param)

    if domain == :spherical
        shape = size(mg.data[1].x)
        msk = get_spherical_domain(shape)

        mskv = get_mask_at_3Dcorners(msk)
        @. mg.oper[1].msk = (param.location == :vertices) ?  mskv : msk

        setup_operators!(mg)
    end


    reset_all_levels(mg)
    b = mg.data[1].b
    msk = mg.oper[1].msk
    nx,ny,nz = size(b)

    if config == :random
        b[:,:,:] .= randn((nx,ny,nz))

    elseif config == :dipole
        di = div(n,4)
        nh = 3
        i0 = nh + div(n,2)
        b[i0-di,i0,i0] = 1
        b[i0+di,i0,i0] = -1

    elseif config == :gradient
        for k in 1:nz, j in 1:ny, i in 1:nx
            b[i,j,k] = (i-3-(nx-6)/2-0.5)*msk[i,j,k]
        end
    else
        @assert false "config not defined"
    end

    @. b *= msk

    subtract_mean!(b, msk)
    res, ite = solve!(mg)

end
