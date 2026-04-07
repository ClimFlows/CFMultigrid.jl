#include("highlevel_operators.jl")


function setup(shape,nhalo,ope,bc;kwargs...)
    levels = setup_levels(shape...;nhalo=nhalo)
    setup_gmg(levels,ope,bc;kwargs...)
end

function setup_levels(nx,ny,nz;nhalo=1)
    maxlevs = 29
    location = (HALOED,HALOED,CLOSED)
    location = (CLOSED,CLOSED,CLOSED)
    axes = Axes(nx,ny,nz,nhalo,location...)
    switch = (nz>1) ? XYZ : XY
    levels = [(axes, switch)]

    lev=1
    while true
        ((nx*ny>4) && (nx>2) && (ny>2)) || break
        (lev<maxlevs) || break
        nx = div(nx,2)
        ny = div(ny,2)
        if switch == XYZ
            nz = div(nz, 2)
        end
        if nz == 1
            switch = XY
        end
        axes = Axes(nx,ny,nz,nhalo,location...)
        push!(levels, (axes, switch))
        lev += 1
    end

    return levels
end

function setup_gmg(levels,ope,bc;kwargs...)
    mg = Tuple([Grid(axes,ope,bc,mode;kwargs...) for (axes, mode) in levels])
    for k in 1:length(mg)-1
        println("LEVEL $k")
        set_Rcoef(mg[k],mg[k+1])
        set_Pcoef(mg[k],mg[k+1])
        set_ope_coef(mg[k],mg[k+1])
    end
    return mg
end

function set_Rcoef(fine, coarse)

    for (i,j,k) in CENTERS(fine.axes)
        fine.b[i,j,k] = 1
    end

    for (i,j,k) in CENTERS(coarse.axes)
        coarse.Rcoef[i,j,k] = 1
    end

    restriction(coarse.b, fine.b, coarse)

    for (i,j,k) in CENTERS(coarse.axes)
        coarse.Rcoef[i,j,k] = coarse.b[i,j,k]>0 ? 4/coarse.b[i,j,k] : 0
    end
    @. fine.b = 0
    @. coarse.b = 0
end

function set_Pcoef(fine,coarse)

    @. fine.x = 0
    @. coarse.x = 0

    for (i,j,k) in CENTERS(fine.axes)
        fine.Pcoef[i,j,k] = 1
    end
    for (i,j,k) in CENTERS(coarse.axes)
        coarse.x[i,j,k] = 1
    end

    prolongation(fine.x, coarse.x, fine, coarse)

    for (i,j,k) in CENTERS(fine.axes)
        fine.Pcoef[i,j,k] = fine.x[i,j,k]>0 ? 1/fine.x[i,j,k] : 0
    end

    @. fine.x = 0
    @. coarse.x = 0

end

