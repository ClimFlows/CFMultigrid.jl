using Printf
#using Plots

include("gmg_types.jl")
include("halo.jl")
include("gmg_setup.jl")
include("gmg_cycles.jl")

function solve!(mg:: Gmg)

    cycles = Dict([ (:onelev, onelevel!),
                    (:twolev, twolevels!),
                    (:nlev, nlevels!),
                    (:v, vcycle!),
                    (:vrec, vcycle_recursive!),
                    (:f, fcycle!)])

    @assert mg.param.cycle in keys(cycles)

    cycle = cycles[mg.param.cycle]

    normb = norm(mg, 1, :b)

    #if mg.param.verbose println("||b||=", normb) end
    res = normresidual!(mg, normb)

    nite = 0
    printstat(mg, nite, res;header=true)
    while true
        ((res>mg.param.tol) && (nite<mg.param.maxite)) || break
        nite += 1
        cycle(mg)
        res = normresidual!(mg, normb)
        # heatmap(mg.data[1].r[:,:,1],show=true,
        #     aspect_ratio=true,colormap=:RdBu)

        # sleep(1.0)
        printstat(mg, nite, res)
    end

    return (res, nite)
end


function printstat(mg:: Gmg, nite, res; header=false)
    if mg.param.verbose
        if header
            println(" ----------------")
            println("  ite     res    ")
            println(" ----------------")
        end
        @printf " %3i   %.3e\n" nite res
    end
end
