#include("grids.jl")

function restriction(xc,xf,coarse::G) where {RP,AX,T,O,G<:Grid{RP,AX,T,O}}
    restriction(xc,xf,coarse.Rcoef,coarse.axes,RP)
end

function prolongation(xf,xc,fine::G,coarse) where {RP,AX,T,O,G<:Grid{RP,AX,T,O}}
    prolongation(xf,xc,fine.Pcoef,coarse.axes,RP)
end

function residual(grid)
    (;r,x,b,ope) = grid
    residual(r,x,b,ope)
    #apply_lateral_bc(r,grid.axes,Neumann)
end


function relativenormresidual(grid::G,normb) where {G<:Grid}
    (;r,x,b,ope,axes) = grid
    residual(r,x,b,ope)
    norm(grid,:r)/normb
end

function norm(grid, which)
    total = 0.
    x = which == :b ? grid.b : grid.r
    (;ax1,ax2,ax3) = grid.axes
    for i in CENTERS(ax1), j in CENTERS(ax2), k in CENTERS(ax3)
        total += x[i,j,k]^2
    end
    return total
end
