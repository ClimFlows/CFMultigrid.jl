function sloop3d!(x:: A, y:: A, b:: A, idiag:: A, grid:: Grid, omega:: Float64)
    for k in kcp(grid), j in jcp(grid), i in icp(grid)
        y[i,j,k] = (1.0-omega)*x[i,j,k]-omega*idiag[i,j,k]*(
            b[i,j,k]-(
                    x[i-1,j,k]+x[i+1,j,k]+
                    x[i,j-1,k]+x[i,j+1,k]+
                    x[i,j,k-1]+x[i,j,k+1]))
    end
end

function sloop3d_v1!(x:: A, y:: A, b:: A, idiag:: A, grid:: Grid, omega:: Float64, n1, n2)
    irange = range(1+n2,length(x)-n2)
    for I in irange
        y[I] = (1.0-omega)*x[I]-omega*idiag[I]*(
            b[I]-(
                    x[I-1]+x[I+1]+
                    x[I-n1]+x[I+n1]+
                    x[I-n2]+x[I+n2]))
    end
end



function smooth_laplacian!(mg:: Gmg, lev:: Int64, nite:: Int64)
    grid = mg.levels[lev].grid

    idiag = mg.oper[lev].idiag
    x =mg.data[lev].x
    y =mg.data[lev].y
    b =mg.data[lev].b

    n1 = size(x,1)
    n2 = n1*size(x,2)

    for kt in 1:nite
        sloop3d_v1!(x,y,b,idiag,grid,mg.param.omega,n1,n2)
        sloop3d_v1!(y,x,b,idiag,grid,mg.param.omega,n1,n2)
        fill!(mg, lev, x)
    end

end

function rloop3d!(x:: A, r:: A, b:: A, diag:: A, msk:: M, grid:: Grid)
    for k in kc(grid), j in jc(grid), i in ic(grid)
        r[i,j,k] = msk[i,j,k]*(diag[i,j,k]*x[i,j,k]+
            b[i,j,k]-(
                    x[i-1,j,k]+x[i+1,j,k]+
                    x[i,j-1,k]+x[i,j+1,k]+
                    x[i,j,k-1]+x[i,j,k+1]))
    end
end

function residual_laplacian!(mg:: Gmg, lev:: Int64)
    grid = mg.levels[lev].grid
    msk = mg.oper[lev].msk
    diag = mg.oper[lev].diag
    x =mg.data[lev].x
    r =mg.data[lev].r
    b =mg.data[lev].b

    rloop3d!(x,r,b,diag,msk,grid)

    fill!(mg, lev, r)
end

@loops function sloop2d!(_, x:: A, y:: A, b:: A, idiag:: A, msk:: M, grid:: Grid, omega:: Float64, n1)
    let Irange = range(1+n1,length(x)-n1)
        @vec for I in Irange
            y[I] = (1.0-omega)*x[I]-omega*idiag[I]*(
                b[I]-(x[I-1]+x[I+1]+x[I-n1]+x[I+n1]))
        end
    end
end


@loops function rloop2d_v1!(_, r:: A, x:: A, b:: A, diag:: A, msk:: M, grid:: Grid, n1)
    let Irange = range(1+n1,length(x)-n1)
        @vec for I in Irange
            r[I] = msk[I]*(diag[I]*x[I]+
                b[I]-(
                    x[I-1]+x[I+1]+
                        x[I-n1]+x[I+n1]))
        end
    end
end

function rloop2d!(r:: A, x:: A, b:: A, diag:: A, msk:: M, grid:: Grid, n1)
        for I in range(1+n1,length(x)-n1)
            r[I] = msk[I]*(diag[I]*x[I]+
                b[I]-(x[I-1]+x[I+1]+x[I-n1]+x[I+n1]))
        end

end

function smooth_laplacian2d!(mg:: Gmg, lev:: Int64, nite:: Int64)
    mgr = mg.mgr[lev]
    grid = mg.levels[lev].grid
    msk = mg.oper[lev].msk
    idiag = mg.oper[lev].idiag
    x =mg.data[lev].x
    y =mg.data[lev].y
    b =mg.data[lev].b

    n1 = size(x,1)

    for kt in 1:nite
        sloop2d!(mgr, x,y,b,idiag,msk,grid,mg.param.omega, n1)
        sloop2d!(mgr, y,x,b,idiag,msk,grid,mg.param.omega, n1)
        fill!(mg, lev, x)
    end

end

function residual_laplacian2d!(mg:: Gmg, lev:: Int64)
    mgr = mg.mgr[lev]
    grid = mg.levels[lev].grid
    msk = mg.oper[lev].msk
    diag = mg.oper[lev].diag
    x =mg.data[lev].x
    r =mg.data[lev].r
    b =mg.data[lev].b

    n1 = size(x,1)

    rloop2d!(r,x,b,diag,msk,grid, n1)
    fill!(mg, lev, r)

end
