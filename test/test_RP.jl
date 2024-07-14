# using AnotherGeometricMultigrid
# using Plots

alloc_like(x) = x .* 0

function test_2d_msk_centers()
    mg, res, ite = test_convergence(128*8,:centers;config=:dipole);
    return isapprox(res, 1.67e-13, atol = 1e-14) & (ite==8)
end

function test_2d_msk_vertices()
    mg, res, ite = test_convergence(128*8,:vertices;config=:dipole);
    return isapprox(res, 1.639e-13, atol = 1e-14) & (ite==8)
end

function test_2d_msk_vertices_convergence_vs_rhs()
    out = true
    for config in [:dipole, :random, :gradient]
        mg, res, ite = test_convergence(256,:vertices;config=config);
        out = out & (ite ==8)
    end
    return out
end

function test_2d_msk_vertices_convergence_vs_resolution()
    out = true
    for p in 7:10
        mg, res, ite = test_convergence(128*p,:vertices;config=:dipole);
        out = out & (ite ==8)
    end
    return out
end

function test_2d_msk_centers_convergence_vs_resolution()
    out = true
    for p in 7:10
        mg, res, ite = test_convergence(128*p,:centers;config=:random);
        out = out & (ite ==8)
    end
    return out
end




function get_mask_at_corners(msk)
    mskv = alloc_like(msk)
    nx, ny,nz = size(msk)
    for k=1:nz, j = 2:ny, i = 2:nx
        fourcells = msk[i, j] + msk[i-1, j] + msk[i, j-1] + msk[i-1, j-1]
        mskv[i, j,k] = (fourcells == 4 ? 1 : 0)
    end
    return mskv
end


function get_closed_domain(shape; nh = 3)
    msk = zeros(UInt8, shape)
    @. msk[nh+1:end-nh, nh+1:end-nh, :] = 1

    return msk
end

function get_circular_domain(shape; nh = 3)
    msk = get_closed_domain(shape; nh = nh)
    nx, ny, nz = shape
    dx = 1.0 / (nx - 2 * nh)
    dy = 1.0 / (ny - 2 * nh)

    for k = 1:nz, j = nh+1:ny-nh+1, i = 1+nh:nx-nh+1
        x = (i - nh - 0.5) * dx
        y = (j - nh - 0.5) * dy
        d2 = (x - 0.5)^2 + (y - 0.5)^2
        r2 = 0.5^2
        msk[i, j,k] = d2 < r2 ? 1 : 0

    end

    return msk
end

function get_banded_domain(shape; nh = 3)
    msk = get_closed_domain(shape; nh = nh)
    nx, ny, nz = shape
    dx = 1.0 / (nx - 2 * nh)
    dy = 1.0 / (ny - 2 * nh)

    for k = 1:nz, j = nh+1:ny-nh+1, i = 1+nh:nx-nh+1
        x = (i - nh - 0.5) * dx
        y = (j - nh - 0.5) * dy
        d2 = (x - 0.5)^2
        r2 = 0.4^2
        msk[i, j,k] = d2 < r2 ? 1 : 0

    end

    return msk
end


function setup_test(;kwargs=Dict())
    param = Param(;kwargs...)
    param.nz=1
    param.case=:twod
    param.cycle=:v
    param.tol=1e-12
    param.ndeepest = 20
    param.npre = 2
    param.npost = 2
    param.omega = 0.92
    param.maxite = 15
    mg = get_gmg(param)

    shape = size(mg.data[1].x)
    msk = get_circular_domain(shape)

    mskv = get_mask_at_corners(msk)
    @. mg.oper[1].msk = (param.location == :vertices) ?  mskv : msk

    setup_operators!(mg)
    return mg
end

function reset_x_on_all_levels(mg)
    for lev in 1:mg.nlevels
        mg.data[lev].x[:,:,:] .= 0
    end
end

function reset_all_levels(mg)
    for lev in 1:mg.nlevels
        mg.data[lev].x[:,:,:] .= 0
        mg.data[lev].y[:,:,:] .= 0
        mg.data[lev].r[:,:,:] .= 0
        mg.data[lev].b[:,:,:] .= 0
    end
end

function mgsetup(;location=:vertices,n=128,verbose=true)
    p = Dict([(:nx,n),
              (:ny,n),
              (:nz,1),
              (:location, location),
              (:verbose, verbose),
              (:case,:twod)])

    mg = setup_test(;kwargs=p)

    reset_all_levels(mg)
    return mg
end

function subtract_mean!(b, msk)
    @. b *= msk

    bmean = sum(b)/sum(msk)

    @. b -= bmean*msk

end

function test_convergence(n,location; config=:dipole, verbose=false)
    mg = mgsetup(;location=location,n=n,verbose=verbose)
    reset_x_on_all_levels(mg)

    b = mg.data[1].b
    msk = mg.oper[1].msk

    nx,ny,nz = size(b)

    if config == :dipole
        i0,j0 = div(nx,2),div(ny,2)
        b[i0-4,j0,1] = -1
        b[i0+4,j0,1] = +1

    elseif config == :gradient
        for k in 1:nz, j in 1:ny, i in 1:nx
            b[i,j,k] = (i-3-(nx-6)/2-0.5)*msk[i,j,k]
        end

    elseif config == :random
        b[:,:,:] .= randn((nx,ny,nz))

    else
        @assert false "config not defined"

    end

    subtract_mean!(b, msk)

    res, ite = solve!(mg)

    if verbose
        out = mg.data[1].x[:,:,:]*1.0
        @. out[msk == 0] = NaN

        heatmap(out[:,:,1],levels=50,show=true,
                aspect_ratio=true,colormap=:RdBu)
    end

    return mg, res, ite
end

function set_b_finest(mg)
    @. mg.data[1].b = 1
    # for I in eachindex(mg.data[1].b)
    #     mg.data[1].b[I] = 1#-mg.oper[1].msk[I]
    # end
    # @. mg.data[1].x = 1
    # mg.residual[1](mg, 1)

end

function set_x_coarsest(mg)
    for lev in 1:mg.nlevels-1
        mg.data[lev].x[:,:,:] .= 0
    end
    mg.data[end].x[:,:,:] .= mg.oper[end].msk
end


function restrict_all(mg)
    for lev in 1:mg.nlevels-1
        mg.restriction[lev](mg, lev, :b)
    end
end

function prolong_all(mg)
    for lev in mg.nlevels-1:-1:1
        mg.prolongation[lev](mg, lev)
    end
end

function test_restrict()
    mg = mgsetup()
    set_b_finest(mg)
    restrict_all(mg)
    heatmap(mg.data[2].b[:,:,1],show=true)
    return mg
end

function test_prolong()
    mg = mgsetup()
    set_x_coarsest(mg)
    prolong_all(mg)
    heatmap(mg.data[1].x[:,:,1],show=true)
    return mg
end
