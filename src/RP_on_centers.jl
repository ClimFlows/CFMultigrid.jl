square(x,i,j,k,di,dj,dk) = return (
    3*(x[i+di,j,k+dk]+x[i,j+dj,k+dk])+
    x[i+di,j+dj,k+dk])

function rcloop3d!(x:: A, xc:: A, coef:: A, grid:: Grid)
    nh = grid.nh

    for kk in kc(grid)
        k = c2fidx(kk, nh)
        for jj in jc(grid)
            j = c2fidx(jj, nh)
            for ii in ic(grid)
                i = c2fidx(ii, nh)

                xc[ii,jj,kk] = coef[ii,jj,kk]*(
                        x[i,j  ,k  ] + x[i+1,j  ,k  ]+
                        x[i,j+1,k  ] + x[i+1,j+1,k  ]+
                        x[i,j  ,k+1] + x[i+1,j  ,k+1]+
                        x[i,j+1,k+1] + x[i+1,j+1,k+1])
            end
        end
    end
end

function global_idx(grid)
    nh = grid.nh
    I0 = nh*grid.tx*grid.ty+nh*grid.ty+nh+1
    I1 = I0+grid.nz*grid.tx*grid.ty
    return I0:I1
end

function fine_idx(Ic, grid)
    nh = grid.nh
    tx = grid.nx*2+2*nh
    ty = grid.ny*2+2*nh

    idx, ii = divrem(Ic-1, grid.tx)
    kk, jj = divrem(idx, grid.ty)

    k = (kk-nh)*2+nh
    j = (jj-nh)*2+nh
    i = (ii-nh)*2+nh

    I = (k*ty+j)*tx+i+1

    return I
end

function rcloop3d_v2!(x:: A, xc:: A, coef:: A, grid:: Grid)
    # implem with a single index loop (anticipating GPU)
    # works but slower, probably because of fine_idx()
    # TODO improve fine_idx()
    nh = grid.nh
    tx = grid.nx*2+2*nh
    ty = grid.ny*2+2*nh

    di, dj, dk = 1, tx, tx*ty

    for Ic in global_idx(grid)
        I = fine_idx(Ic, grid)

        xc[Ic] = coef[Ic]*(
            x[I] +
                x[I+di] + x[I+dj] + x[I+dk] +
                x[I+di+dj] + x[I+di+dk] + x[I+dj+dk] +
                x[I+di+dj+dk] )
    end
end




function restriction_centers!(mg:: Gmg, lev:: Int64, field:: Symbol)
    grid = mg.levels[lev+1].grid
    x = getfield(mg.data[lev], field)
    xc = mg.data[lev+1].b
    coef = mg.oper[lev+1].Rcoef

    rcloop3d!(x, xc, coef, grid)

end

function pcloop3d!(x:: A, xf:: A, coef:: A, grid:: Grid)
    nh = grid.nh
    for k in nh+1:2:grid.nz+nh
        kc = f2cidx(k, nh)
        for j in nh+1:2:grid.ny+nh
            jc = f2cidx(j, nh)
            for i in nh+1:2:grid.nx+nh
                ic = f2cidx(i, nh)
                x28 = x[ic,jc,kc]*28
                @inline for dj in 0:1, di in 0:1
                    a = square(x,ic,jc,kc,2*di-1,2*dj-1,-1)
                    b = square(x,ic,jc,kc,2*di-1,2*dj-1,+0)*3
                    c = square(x,ic,jc,kc,2*di-1,2*dj-1,+1)
                    xf[i+di,j+dj,k  ] += coef[i+di,j+dj,k  ]*(b+a+x28)
                    xf[i+di,j+dj,k+1] += coef[i+di,j+dj,k+1]*(b+c+x28)
                end
            end
        end
    end

end

function prolongation_centers!(mg:: Gmg, lev:: Int64)
    grid = mg.levels[lev].grid
    coef = mg.oper[lev].Pcoef
    xf = mg.data[lev].x
    x = mg.data[lev+1].x

    pcloop3d!(x, xf, coef, grid)
end


function rcloop2d!(x:: A, xc:: A, coef:: A, grid:: Grid)
    nh = grid.nh

    for kk in kc(grid)
        k = kk
        for jj in jc(grid)
            j = c2fidx(jj, nh)
            for ii in ic(grid)
                i = c2fidx(ii, nh)

                xc[ii,jj,kk] = coef[ii,jj,kk]*(
                        x[i,j  ,k  ] + x[i+1,j  ,k  ]+
                        x[i,j+1,k  ] + x[i+1,j+1,k  ])
            end
        end
    end
end

function restriction_centers2d!(mg:: Gmg, lev:: Int64, field:: Symbol)
    grid = mg.levels[lev+1].grid
    x = getfield(mg.data[lev], field)
    xc = mg.data[lev+1].b
    coef = mg.oper[lev+1].Rcoef

    rcloop2d!(x, xc, coef, grid)

end

function pcloop2d!(x:: A, xf:: A, coef:: A, grid:: Grid)
    nh = grid.nh

    for k in kc(grid)
        kc = k
        for j in nh+1:2:grid.ny+nh
            jc = f2cidx(j, nh)
            for i in nh+1:2:grid.nx+nh
                ic = f2cidx(i, nh)
                x9 = x[ic,jc,kc]*9
                @inline for dj in 0:1, di in 0:1
                    b = square(x,ic,jc,kc,2*di-1,2*dj-1,0)
                    xf[i+di,j+dj,k  ] += coef[i+di,j+dj,k  ]*(b+x9)
                end
            end
        end
    end
end

function prolongation_centers2d!(mg:: Gmg, lev:: Int64)
    grid = mg.levels[lev].grid
    coef = mg.oper[lev].Pcoef

    xf = mg.data[lev].x
    x = mg.data[lev+1].x

    pcloop2d!(x, xf, coef, grid)

end
