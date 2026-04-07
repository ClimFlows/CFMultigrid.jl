using Printf

abstract type CYCLE end
abstract type ONELEVEL<:CYCLE end
abstract type TWOLEVEL<:CYCLE end
abstract type VCYCLE<:CYCLE end
abstract type FCYCLE<:CYCLE end

include("setup.jl")

struct Param{I,F,S}
    npre::I
    npost::I
    ndeepest::I
    maxite::I
    tol::F
    smoother::S
end

function solve!(mg, param, ::Type{C};verbose=true) where {C<:CYCLE}
    (;maxite, tol) = param

    finest = mg[1]
    nite = 0

    normb = norm(finest, :b)

    if normb == 0
        fill!(finest.x,0)
        res = 0.0
        return (res, nite)
    end

    res = relativenormresidual(finest, normb)


    log(nite,res,verbose)
    while (res>tol) & (nite<maxite)
        nite += 1
        vcycle(mg, param)
        res = relativenormresidual(finest, normb)
        log(nite,res,verbose)
    end

    return (res, nite)

end

function vcycle(mg::MG, param; k0=1) where {MG<:Tuple{Vararg{Grid}}}
    (;npre,npost,ndeepest,smoother) = param
    levels = k0:length(mg)-1
    for k in levels
        smooth(npre, mg[k], smoother)
        residual(mg[k])
        restriction(mg[k+1].b, mg[k].r, mg[k+1])
        fill!(mg[k+1].x,0)
    end

    smooth(ndeepest, mg[end], smoother)

    for k in reverse(levels)
        prolongation(mg[k].x, mg[k+1].x, mg[k], mg[k+1])
        smooth(npost, mg[k], smoother)
    end

end

log(nite,res,verbose) = verbose ? (@show nite, res) : nothing


Base.show(io::IO, rp::R) where {D,B,R<:RP{D,B}}= print(io,"$D")

Base.show(io::IO, grid::G) where {G<:Grid} = begin
    (;axes)=grid
    (;ax1,ax2,ax3) = axes
    f = Printf.Format("%4i x %4i x %4i")
    Printf.format(io, f, ax1.n, ax2.n, ax3.n)
end

Base.show(io::IO, mg::MG) where {MG<:Tuple{Vararg{Grid}}} = begin
    println(io, "Multigrid:")
    grid = mg[1]
    bc = typeof(grid).parameters[1].parameters[2]
    println(io, "  Operator: ", typeof(grid.ope).name.name)
    println(io, "  Boundary conditions: ", bc)

    for (k,grid) in enumerate(mg)
        print(io,"    - ")
        show(grid)
        if k<length(mg)
            D = typeof(grid).parameters[1].parameters[1]
            print(io, "  ",D)
        end
        print(io,"\n")
    end
end
