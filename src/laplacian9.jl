
function smooth_laplacian9!(mg:: Gmg, lev:: Int64, nite:: Int64)
    grid = mg.levels[lev].grid

    idiag = mg.oper[lev].idiag
    x =mg.data[lev].x
    y =mg.data[lev].y
    b =mg.data[lev].b

    omega = mg.param.omega

    cff = 1.0 - omega

    for k in kcp(grid), j in jcp(grid), i in icp(grid)
        y[i,j,k] = cff*x[i,j,k]-omega*idiag[i,j,k]*(
            b[i,j,k]
            -(
                    x[i-1,j,k]+x[i+1,j,k]+
                    x[i,j-1,k]+x[i,j+1,k]+
                    x[i,j,k-1]+x[i,j,k+1] ) / 2
            -(
                x[i-1,j-1,k]+x[i+1,j-1,k]+
                x[i-1,j+1,k]+x[i+1,j+1,k]+
                x[i-1,j,k-1]+x[i+1,j,k-1]+
                x[i-1,j,k+1]+x[i+1,j,k+1]+
                x[i,j-1,k-1]+x[i,j+1,k-1]+
                x[i,j-1,k+1]+x[i,j+1,k+1]
            )/4)
    end

    for k in kc(grid), j in jc(grid), i in ic(grid)
        x[i,j,k] = cff*y[i,j,k]-omega*idiag[i,j,k]*(
            b[i,j,k]-(
                    y[i-1,j,k]+y[i+1,j,k]+
                    y[i,j-1,k]+y[i,j+1,k]+
                    y[i,j,k-1]+y[i,j,k+1])/2
            -(
                y[i-1,j-1,k]+y[i+1,j-1,k]+
                y[i-1,j+1,k]+y[i+1,j+1,k]+
                y[i-1,j,k-1]+y[i+1,j,k-1]+
                y[i-1,j,k+1]+y[i+1,j,k+1]+
                y[i,j-1,k-1]+y[i,j+1,k-1]+
                y[i,j-1,k+1]+y[i,j+1,k+1]
            )/4)
    end

end

function residual_laplacian9!(mg:: Gmg, lev:: Int64)
    grid = mg.levels[lev].grid
    msk = mg.oper[lev].msk
    diag = mg.oper[lev].diag
    x =mg.data[lev].x
    r =mg.data[lev].r
    b =mg.data[lev].b

    for k in kc(grid), j in jc(grid), i in ic(grid)
        r[i,j,k] = msk[i,j,k]*(diag[i,j,k]*x[i,j,k]+
            b[i,j,k]-(
                    x[i-1,j,k]+x[i+1,j,k]+
                    x[i,j-1,k]+x[i,j+1,k]+
                    x[i,j,k-1]+x[i,j,k+1])/2
            -(
                x[i-1,j-1,k]+x[i+1,j-1,k]+
                x[i-1,j+1,k]+x[i+1,j+1,k]+
                x[i-1,j,k-1]+x[i+1,j,k-1]+
                x[i-1,j,k+1]+x[i+1,j,k+1]+
                x[i,j-1,k-1]+x[i,j+1,k-1]+
                x[i,j-1,k+1]+x[i,j+1,k+1]
            )/4)

    end
end
