include("operators_at_level.jl")

function onelevel!(mg:: Gmg)
    lev = 1
    mg.smoother[lev](mg, lev, mg.param.ndeepest)
end

function twolevels!(mg:: Gmg)

    lev = 1
    mg.smoother[lev](mg, lev, mg.param.npre)
    mg.residual[lev](mg, lev)
    mg.restriction[lev](mg, lev, :r)
    set_to_zero!(mg, lev+1, :x)

    lev = 2
    mg.smoother[lev](mg, lev, mg.param.ndeepest)

    lev = 1
    mg.prolongation[lev](mg, lev)
    mg.smoother[lev](mg, lev, mg.param.npost)

end

function nlevels!(mg:: Gmg; nlev=3)

    for lev in 1:nlev-1
        mg.smoother[lev](mg, lev, mg.param.npre)
        mg.residual[lev](mg, lev)
        mg.restriction[lev](mg, lev, :r)
        set_to_zero!(mg, lev+1, :x)
    end

    mg.smoother[nlev](mg, mg.nlevels, mg.param.ndeepest)

    for lev in nlev-1:1
        mg.prolongation[lev](mg, lev)
        mg.smoother[lev](mg, lev, mg.param.npost)
    end

end

function vcycle_recursive!(mg:: Gmg; lev=1)
    if lev == mg.nlevels
        mg.smoother[mg.nlevels](mg, mg.nlevels, mg.param.ndeepest)
    else
        mg.smoother[lev](mg, lev, mg.param.npre)
        mg.residual[lev](mg, lev)
        if mg.param.verbose
            res1 = norm(mg, lev, :r)
        end
        mg.restriction[lev](mg, lev, :r)
        set_to_zero!(mg, lev+1, :x)

        vcycle_recursive!(mg;lev=lev+1)

        mg.prolongation[lev](mg, lev)
        mg.smoother[lev](mg, lev, mg.param.npost)
        if mg.param.verbose
            mg.residual[lev](mg, lev)
            res2 = norm(mg, lev, :r)
            @info lev, res1/res2
        end

    end
end
function vcycle!(mg:: Gmg; lev0=1)

    for lev in lev0:mg.nlevels-1
        mg.smoother[lev](mg, lev, mg.param.npre)
        mg.residual[lev](mg, lev)
        mg.restriction[lev](mg, lev, :r)
        set_to_zero!(mg, lev+1, :x)
    end

    mg.smoother[mg.nlevels](mg, mg.nlevels, mg.param.ndeepest)

    for lev in mg.nlevels-1:-1:lev0
        mg.prolongation[lev](mg, lev)
        mg.smoother[lev](mg, lev, mg.param.npost)
    end

end

function fcycle!(mg:: Gmg)

    for lev in 1:mg.nlevels-1
        field = (lev>1) ? :b : :r
        mg.restriction[lev](mg, lev, field)
        set_to_zero!(mg, lev+1, :x)
    end

    mg.smoother[mg.nlevels](mg, mg.nlevels, mg.param.ndeepest)

    for lev in mg.nlevels-1:-1:1
        mg.prolongation[lev](mg, lev)
        vcycle!(mg; lev0=lev)
    end

end
