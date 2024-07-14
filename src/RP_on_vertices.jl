function restriction_vertices!(mg:: Gmg, lev:: Int64, field:: Symbol)
    grid = mg.levels[lev+1].grid
    x = getfield(mg.data[lev], field)
    xc = mg.data[lev+1].b
    coef = mg.oper[lev+1].Rcoef
    nh = grid.nh

    for kk in kc(grid)
        k = c2fidx(kk, nh)
        for jj in jc(grid)
            j = c2fidx(jj, nh)
            for ii in ic(grid)
                i = c2fidx(ii, nh)

                xc[ii,jj,kk] = coef[ii,jj,kk]*(
                  +  x[i-1,j-1,k-1]+2*x[i,j-1,k-1]+  x[i+1,j-1,k-1]
                  +2*x[i-1,j  ,k-1]+4*x[i,j  ,k-1]+2*x[i+1,j  ,k-1]
                  +  x[i-1,j+1,k-1]+2*x[i,j+1,k-1]+  x[i+1,j+1,k-1]
                    #
                  +2*x[i-1,j-1,k]+4*x[i,j-1,k]+2*x[i+1,j-1,k]
                  +4*x[i-1,j  ,k]+8*x[i,j  ,k]+4*x[i+1,j  ,k]
                  +2*x[i-1,j+1,k]+4*x[i,j+1,k]+2*x[i+1,j+1,k]
                    #
                  +  x[i-1,j-1,k+1]+2*x[i,j-1,k+1]+  x[i+1,j-1,k+1]
                  +2*x[i-1,j  ,k+1]+4*x[i,j  ,k+1]+2*x[i+1,j  ,k+1]
                  +  x[i-1,j+1,k+1]+2*x[i,j+1,k+1]+  x[i+1,j+1,k+1]
                )
            end
        end
    end
end

function prolongation_vertices!(mg:: Gmg, lev:: Int64)
    grid = mg.levels[lev].grid
    coef = mg.oper[lev].Pcoef
    nh = grid.nh

    xf = mg.data[lev].x
    x = mg.data[lev+1].x

    for k in nh+1:2:grid.nz+nh
        kc = f2cidx(k, nh)
        for j in nh+1:2:grid.ny+nh
            jc = f2cidx(j, nh)
            for i in nh+1:2:grid.nx+nh
                ic = f2cidx(i, nh)

                xf[i  ,j  ,k  ] += coef[i  ,j  ,k  ]*x[ic,jc,kc]
                xf[i+1,j  ,k  ] += coef[i+1,j  ,k  ]*(x[ic,jc,kc]+x[ic+1,jc,kc])
                xf[i  ,j+1,k  ] += coef[i  ,j+1,k  ]*(x[ic,jc,kc]+x[ic,jc+1,kc])
                xf[i  ,j  ,k+1] += coef[i  ,j  ,k+1]*(x[ic,jc,kc]+x[ic,jc,kc+1])

                xf[i+1,j+1,k  ] += coef[i+1,j+1,k  ]*(
                    x[ic,jc,kc]+x[ic+1,jc,kc]+x[ic,jc+1,kc]+x[ic+1,jc+1,kc])

                xf[i+1,j  ,k+1] += coef[i+1,j  ,k+1]*(
                    x[ic,jc,kc]+x[ic+1,jc,kc]+x[ic,jc,kc+1]+x[ic+1,jc,kc+1])

                xf[i  ,j+1,k+1] += coef[i  ,j+1,k+1]*(
                    x[ic,jc,kc]+x[ic,jc,kc+1]+x[ic,jc+1,kc]+x[ic,jc+1,kc+1])

                xf[i+1,j+1,k+1] += coef[i+1,j+1,k+1]*(
                   x[ic,jc,kc]+x[ic+1,jc,kc]+x[ic,jc,kc+1]+x[ic+1,jc,kc+1]+
                    x[ic,jc+1,kc]+x[ic+1,jc+1,kc]+x[ic,jc+1,kc+1]+x[ic+1,jc+1,kc+1])

            end
        end
    end

end

function rvloop2d_vorig!(x:: A, xc:: A, coef:: A, grid:: Grid)
    nh = grid.nh

    for kk in kc(grid)
        k = kk
        for jj in nh+1:grid.ny+nh#jc(grid)#1+grid.nh+1:grid.ny+grid.nh-1#jv(grid)
            j = c2fidx(jj, nh)
            for ii in nh+1:grid.nx+nh#ic(grid)# 1+grid.nh+1:grid.nx+grid.nh-1#iv(grid)
                i = c2fidx(ii, nh)

                xc[ii,jj,kk] = coef[ii,jj,kk]*(
                  +  x[i-1,j-1,k]+2*x[i,j-1,k]+  x[i+1,j-1,k]
                  +2*x[i-1,j  ,k]+4*x[i,j  ,k]+2*x[i+1,j  ,k]
                  +  x[i-1,j+1,k]+2*x[i,j+1,k]+  x[i+1,j+1,k]
                )

                # xc[ii,jj,kk] = coef[ii,jj,kk]*x[i,j  ,k]*16
            end
        end
    end
end

function rvloop2d!(x:: A, xc:: A, coef:: A, grid:: Grid)
    nh = grid.nh

    for kk in kc(grid)
        k = kk
        for j in nh+1:2:2*grid.ny+nh
            jc = f2cidx(j, nh)
            for i in nh+1:2:2*grid.nx+nh
                ic = f2cidx(i, nh)

                xc[ic,jc,kk] = coef[ic,jc,kk]*(
                  +  x[i-1,j-1,k]+2*x[i,j-1,k]+  x[i+1,j-1,k]
                  +2*x[i-1,j  ,k]+4*x[i,j  ,k]+2*x[i+1,j  ,k]
                  +  x[i-1,j+1,k]+2*x[i,j+1,k]+  x[i+1,j+1,k]
                )

            end
        end
    end
end

function pvloop2d!(x:: A, xf:: A, coef:: A, grid:: Grid)
    nh = grid.nh

    for k in kc(grid)
        kc = k
        for j in nh+1:2:grid.ny+nh
            jc = f2cidx(j, nh)
            for i in nh+1:2:grid.nx+nh
                ic = f2cidx(i, nh)

                xf[i,  j  ,k] += coef[i  ,j  ,k]*(x[ic,jc,kc])
                xf[i+1,j  ,k] += coef[i+1,j  ,k]*(x[ic,jc,kc]+x[ic+1,jc,kc])
                xf[i  ,j+1,k] += coef[i  ,j+1,k]*(x[ic,jc,kc]+x[ic,jc+1,kc])

                xf[i+1,j+1,k] += coef[i+1,j+1,k]*(x[ic,jc,kc]+x[ic+1,jc,kc]
                                                 +x[ic,jc+1,kc]+x[ic+1,jc+1,kc])

            end
        end
    end
end


function restriction_vertices2d!(mg:: Gmg, lev:: Int64, field:: Symbol)
    grid = mg.levels[lev+1].grid
    x = getfield(mg.data[lev], field)
    xc = mg.data[lev+1].b
    coef = mg.oper[lev+1].Rcoef
    nh = grid.nh

    rvloop2d!(x, xc, coef, grid)

end

function prolongation_vertices2d!(mg:: Gmg, lev:: Int64)
    grid = mg.levels[lev].grid
    coef = mg.oper[lev].Pcoef
    nh = grid.nh

    xf = mg.data[lev].x
    x = mg.data[lev+1].x

    pvloop2d!(x, xf, coef, grid)
end
