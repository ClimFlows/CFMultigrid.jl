using ManagedLoops: @loops, @vec
include("laplacian.jl")
include("laplacian9.jl")
include("RP_on_centers.jl")
include("RP_on_vertices.jl")

function set_to_zero!(mg:: Gmg, lev, field:: Symbol)
    x = getfield(mg.data[lev], field)
    x[:,:,:] .= 0.0
end

function norm(mg:: Gmg, lev:: Int64, field:: Symbol)
    grid = mg.levels[lev].grid
    msk = mg.oper[lev].msk
    x = getfield(mg.data[lev], field)

    sum = 0.0
    for k in ka(grid), j in ja(grid), i in ia(grid)
        sum += msk[i,j,k]*x[i,j,k]^2
    end

    if sum == 0.0
        sum = 1.0
    end

    return sum
end

function normresidual!(mg:: Gmg, normb:: Float64)
    if normb > 0.0
        mg.residual[1](mg, 1)
        return norm(mg, 1, :r)/normb
    else
        return 0.0
    end
end


