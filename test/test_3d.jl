include("test_RP.jl")

function test_3d_centers()
    res, ite = run_3d_test(128,:centers;domain=:square,config=:dipole)
    return isapprox(res, 2.55e-13, atol = 1e-14) & (ite==6)
end

function test_3d_vertices()
    res, ite = run_3d_test(128,:vertices;domain=:spherical,config=:dipole)
    return isapprox(res, 2.60e-13, atol = 1e-14) & (ite==6)
end



function get_mask_at_3Dcorners(msk)
    mskv = alloc_like(msk)
    nx, ny,nz = size(msk)
    for k=2:nz, j = 2:ny, i = 2:nx
        eightcells = msk[i, j,k] + msk[i-1, j,k] + msk[i, j-1,k] + msk[i-1, j-1,k]+msk[i, j,k-1] + msk[i-1, j,k-1] + msk[i, j-1,k-1] + msk[i-1, j-1,k-1]
        mskv[i, j,k] = (eightcells == 8 ? 1 : 0)
    end
    return mskv
end


function get_closed_3Ddomain(shape; nh = 3)
    msk = zeros(UInt8, shape)
    @. msk[nh+1:end-nh, nh+1:end-nh, nh+1:end-nh] = 1

    return msk
end

function get_spherical_domain(shape; nh = 3)
    msk = get_closed_domain(shape; nh = nh)
    nx, ny, nz = shape
    dx = 1.0 / (nx - 2 * nh)
    dy = 1.0 / (ny - 2 * nh)
    dz = 1.0 / (nz - 2 * nh)

    for k = nh+1:nz-nh+1, j = nh+1:ny-nh+1, i = 1+nh:nx-nh+1
        x = (i - nh - 0.5) * dx
        y = (j - nh - 0.5) * dy
        z = (k - nh - 0.5) * dz
        d2 = (x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2
        r2 = 0.5^2
        msk[i,j,k] = d2 < r2 ? 1 : 0

    end

    return msk
end

function setup_3d(domain;kwargs=Dict())
    param = Param(;kwargs...)
    param.nz=param.nx
    param.case=:threed
    param.cycle=:v
    param.tol=1e-12
    param.ndeepest = 20
    param.npre = 2
    param.npost = 2
    param.omega = 0.88
    param.omega = 0.84
    param.maxite = 15
    mg = get_gmg(param)

    if domain == :spherical
        shape = size(mg.data[1].x)
        msk = get_spherical_domain(shape)

        mskv = get_mask_at_3Dcorners(msk)
        @. mg.oper[1].msk = (param.location == :vertices) ?  mskv : msk

        setup_operators!(mg)
    end

    return mg
end

function run_3d_test(n,location;
                     domain = :square,config=:dipole,verbose=false)
    p = Dict([(:nx,n),
              (:ny,n),
              (:nz,n),
              (:location, location),
              (:verbose, verbose),
              (:case,:threed)])
    mg = setup_3d(domain;kwargs=p)
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

    # heatmap(mg.data[1].r[:,:,i0] ,show=true)
    # heatmap(mg.data[1].r[:,i0,:]*1e6 ,show=true)
end
