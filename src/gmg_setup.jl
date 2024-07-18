using Printf
using LoopManagers:
    SIMD, PlainCPU, VectorizedCPU, MultiThread, KernelAbstractions_GPU as GPU
# main user function

"""
    mg = get_gmg(param)
Returns a multigrid object, whose specs are defined by param (see ?Param)
  - The operator A is defined by param.operator. So far only the `:laplacian` is implemented
  - Two types of boundary conditions (BC) are possible: Dirichlet and Neumann. They are selected by setting the param.location of the discretized variables: either `:vertices` (for Dirichlet) or `:centers` (for Neumann).

To see the grid hierarchy

    `print(mg.levels)`

To solve A * x = b

    `res, nite = solve!(mg)`

   - b, the r.h.s., must be copied in `mg.data[1].b` (a 3d array)
   - x, the first guess and the solution, is in `mg.data[1].x` (a 3d array)

Irregular domains are handled with a mask system. Look inside 'tests/' for examples.

"""
function get_gmg(param:: Param)
    levels = create_levels(param)
    mg = allocate_mg(param, levels)
    set_default_msk!(mg)
    setup_operators!(mg)
    for lev in 1:mg.nlevels, var in [:b, :x, :y, :r]
        set_to_zero!(mg, lev, var)
    end
    return mg
end


# internal functions

function create_levels(param:: Param)
    nx = param.nx
    ny = param.ny
    nz = param.nz
    nh = param.nh
    case = param.case
    maxlevs = 29#param.maxlevs
    grid = Grid(nx,ny,nz,nh,case)
    levels = [Level(grid)]

    lev=1
    while true
        ((nx>4) && (ny>4)) || break
        (lev<maxlevs) || break
        nx = div(nx,2)
        ny = div(ny,2)
        if case == :threed
            nz = div(nz, 2)
        end
        grid = Grid(nx,ny,nz,nh,case)
        push!(levels, Level(grid))
        lev += 1
    end
    levels
end

Base.show(io::IO, levels::Vector{Level}) = begin
    println("---------------------")
    println(" lev  nx    ny    nz ")
    println("---------------------")
    for (k, level) in enumerate(levels)
        g = level.grid
        @printf " %2i  %3i x %3i x %3i\n" k g.nx g.ny g.nz
    end
end


function allocate_mg(param:: Param, levels:: Vector{Level})
    nlevels = length(levels)
    data = [Data(levels[i].grid) for i in 1:nlevels]
    oper = [Operators(levels[i].grid) for i in 1:nlevels]
    smoothers = []
    residuals = []
    restrictions = []
    prolongations = []
    managers = []
    for lev in 1:nlevels
        smo, res = get_smoother_and_residual(param, levels[lev].grid)
        restri, prolong = get_restriction_and_prolongation(param, levels[lev].grid)
        # if lev<=2
        #     mgr = MultiThread(VectorizedCPU(8))
        # else
        #     mgr = PlainCPU()
        # end
        mgr = PlainCPU()

        push!(smoothers, smo)
        push!(residuals, res)
        push!(restrictions, restri)
        push!(prolongations, prolong)
        push!(managers, mgr)
    end
    Gmg(managers, nlevels, param, levels, data, oper,
        smoothers, residuals, restrictions, prolongations)
end

function set_default_msk!(mg:: Gmg;lev=1)
    vertices = mg.param.location == :vertices
    msk = mg.oper[1].msk
    grid = mg.levels[1].grid
    nh = grid.nh
    if vertices
        kidx = (mg.param.case == :threed) ? (nh+1:nh+grid.nz) : (1:grid.nz)
        for k in kidx, j in jv(grid), i in iv(grid)
            msk[i,j,k] = 1
        end
    else
        for k in kc(grid), j in jc(grid), i in ic(grid)
            msk[i,j,k] = 1
        end
    end
    fill!(mg,lev,msk)
end

function get_smoother_and_residual(param:: Param, grid:: Grid)
    if param.case == :threed
        if param.operator == :laplacian
            return (smooth_laplacian!, residual_laplacian!)
        elseif param.operator == :laplacian9
            return (smooth_laplacian9!, residual_laplacian9!)
        end
        @assert false
    else
        if param.operator == :laplacian
            return (smooth_laplacian2d!, residual_laplacian2d!)
        end
        @assert false
    end

end

function get_restriction_and_prolongation(param:: Param, grid:: Grid)
    if param.case == :threed
        if param.location == :centers
            return (restriction_centers!, prolongation_centers!)
        elseif param.location == :vertices
            return (restriction_vertices!, prolongation_vertices!)
        end
        @assert false
    else
        if param.location == :centers
            return (restriction_centers2d!, prolongation_centers2d!)
        elseif param.location == :vertices
            return (restriction_vertices2d!, prolongation_vertices2d!)
        end
        @assert false
    end
end

function compute_msk!(mg:: Gmg, lev:: Int64)
    vertices = mg.param.location == :vertices
    grid = mg.levels[lev+1].grid
    msk = mg.oper[lev+1].msk
    Rcoef = mg.oper[lev+1].Rcoef
    b = mg.data[lev+1].b

    #coef = (mg.param.case == :threed) ? 2.0 : 2.0
    threshold = vertices ? 10 : 0

    msk[:,:,:] .= 1
    Rcoef[:,:,:] .= 1.0
    mg.data[lev].b[:,:,:] .= mg.oper[lev].msk[:,:,:]

    mg.restriction[lev](mg, lev, :b)
    msk_def = get_default_msk!(mg:: Gmg;lev=lev+1)
    @. b *= msk_def
    # heatmap(b[:,:,1],show=true)
    # sleep(2)
    # @show unique(b)
    for k in ka(grid), j in ja(grid), i in ia(grid)
        msk[i,j,k] = (b[i,j,k]>threshold) ? 1 : 0
    end

     apply_default_msk!(mg, lev)
end

function compute_Rcoef!(mg:: Gmg, lev:: Int64)
    vertices = mg.param.location == :vertices
    grid = mg.levels[lev+1].grid

    Rcoef = mg.oper[lev+1].Rcoef
    msk = mg.oper[lev+1].msk
    b = mg.data[lev+1].b

    Rcoef .= msk

    mg.data[lev].b .= (vertices ? mg.oper[lev].msk : 1.0)
    #mg.data[lev].b .= 1
    mg.restriction[lev](mg, lev, :b)

    # TODO CHECK for the 3D case
    coef = (mg.param.case == :threed) ? 1/16 : 1/4

    fvertices(b) = begin (b>0.0) ? coef : 0.0 end
    fcenters(b) = begin (b>0.0) ? (4/b) : 0.0 end

    fRcoef = (vertices ? fvertices : fcenters)

    for k in ka(grid), j in ja(grid), i in ia(grid)
        Rcoef[i,j,k] = msk[i,j,k]*fRcoef(b[i,j,k])
    end

end

function compute_Pcoef!(mg:: Gmg, lev:: Int64)
    vertices = mg.param.location == :vertices
    grid = mg.levels[lev].grid

    Pcoef = mg.oper[lev].Pcoef
    msk = mg.oper[lev].msk
    mskc = mg.oper[lev+1].msk
    x = mg.data[lev].x
    xc = mg.data[lev+1].x

    Pcoef .= msk
    x .= 0.0
    xc .= vertices ? 1.0 : mskc
    #xc .= mskc

    mg.prolongation[lev](mg, lev)

    fP(x) = begin (x>0) ? 1/x : 0.0 end

    for k in ka(grid), j in ja(grid), i in ia(grid)
        Pcoef[i,j,k] = fP(x[i,j,k])*msk[i,j,k]
    end

end

function compute_diag!(mg:: Gmg, lev:: Int64)
    vertices = mg.param.location == :vertices

    grid = mg.levels[lev].grid
    x = mg.data[lev].x
    b = mg.data[lev].b
    r = mg.data[lev].r

    diag = mg.oper[lev].diag
    idiag = mg.oper[lev].idiag
    msk = mg.oper[lev].msk

    x .= vertices ? 1.0 : msk
    b .= 0.0
    diag .= 0.0

    mg.residual[lev](mg, lev)

    @. diag = -r*msk
    for k in ka(grid), j in ja(grid), i in ia(grid)
        idiag[i,j,k] = (diag[i,j,k]>0.0) ? 1.0/diag[i,j,k] : 0.0
    end

end

function setup_operators!(mg)
    if false
        # old algo for setting up coarser masks
        # => makes the convergence depends on the resolution
        # residual stalls because of the boundary
        # can be overcome by increasing relaxation passes
        # but this is hiding the dust under the carpet
        for lev in 1:mg.nlevels-1
            compute_msk!(mg, lev)
        end
    else
        # new algo: fixes previous defaults
        # TODO: adapt it to the 3D case
        compute_msk_dev(mg)
    end

    for lev in 1:mg.nlevels-1
        compute_Rcoef!(mg, lev)
        compute_Pcoef!(mg, lev)
    end
    for lev in 1:mg.nlevels
        compute_diag!(mg, lev)
    end
end

function get_default_msk!(mg:: Gmg;lev=1::Int64)
    vertices = mg.param.location == :vertices
    grid = mg.levels[lev].grid
    msk = zeros(UInt8, (grid.tx, grid.ty, grid.tz))
    nh = grid.nh
    nx = grid.nx
    ny = grid.ny
    if vertices
        kidx = (mg.param.case == :threed) ? (nh+1:nh+grid.nz) : (1:grid.nz)
        for k in kidx, j in jv(grid), i in iv(grid)
            msk[i,j,k] = 1
        end
    else
        for k in kc(grid), j in jc(grid), i in ic(grid)
            msk[i,j,k] = 1
        end
    end
    fill!(mg, lev, msk)
    return msk
end

function apply_default_msk!(mg:: Gmg, lev:: Int64)
    vertices = mg.param.location == :vertices
    msk = mg.oper[lev].msk
    #mskcopy = deepcopy(msk)
    msk_def = get_default_msk!(mg:: Gmg;lev=lev)
    #@show msk, msk_def
    @. msk *= msk_def
    #msk .= mskcopy .* msk

end


function compute_msk_dev(mg)
    vertices = mg.param.location == :vertices

    for lev in 1:mg.nlevels-1

        grid = mg.levels[lev+1].grid
        msk = mg.oper[lev+1].msk
        Rcoef = mg.oper[lev+1].Rcoef
        b = mg.data[lev+1].b

        #coef = (mg.param.case == :threed) ? 2.0 : 2.0
        threshold = vertices ? 0.67 : 0.18

        msk0 = copy(msk)
        msk[:,:,:] .= 1
        Rcoef[:,:,:] .= 1.0
        if lev == 1
            mg.data[lev].b[:,:,:] .= mg.oper[lev].msk[:,:,:]
        end

        mg.restriction[lev](mg, lev, :b)

        msk[:] .= msk0[:]

        msk_def = get_default_msk!(mg;lev=lev+1)
        if mg.param.case == :threed
            coef = vertices ? 1/64 : 1/8
        else
            coef = vertices ? 1/16 : 1/4
        end

        @. b *= msk_def*coef

        # k0 = (mg.param.case == :threed) ? div(grid.nz,2)+grid.nh : 1
        # heatmap(b[:,:,k0],show=true)
        # sleep(2)
        # @show unique(b)

        for k in ka(grid), j in ja(grid), i in ia(grid)
            msk[i,j,k] = (b[i,j,k]>threshold) ? 1 : 0
        end

        apply_default_msk!(mg, lev)
    end
end
