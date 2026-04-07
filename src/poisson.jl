struct Poisson{T,AX,M,D} <: OPERATOR
    axes:: AX
    diag:: Array{T,3}
    idiag:: Array{T,3}
    msk::Array{Int8,3}
end

function Poisson(axes::AX,::Type{B}, ::Type{D};msk=nothing) where {AX<:AXES, D<:DIRECTIONS, B<:BC}
    T = Float64
    shape = size(axes,(CENTERS,CENTERS,CENTERS))
    diag, idiag = [zeros(T,shape) for _ in 1:2]
    m, status = check_msk(msk, shape)
    ope = Poisson{T,AX,status,D}(axes,diag,idiag,m)
    set_poisson_diag(ope,B)
    ope
end

function set_poisson_diag(ope::Poisson{T,AX,M,D}, ::Type{B}) where {T,AX,M,D<:DIRECTIONS,B<:BC}
    (;diag,idiag,msk) = ope
    shape = size(diag)
    x,r,b = [zeros(T,shape) for _ in 1:3]
    @. x = 1
    coef = diagcoef(shape)
    @. diag = coef
    residual(r,x,b,ope)

    if B==NEUMANN
        @. diag -= r
    else
        @. diag += r
    end

    if M==ON
        @. idiag[msk==1]  =  1/diag[msk==1]
    else
        @. idiag  =  1/diag
    end
end



function residual(r::A,x::A,b::A,ope::O) where {T,AX,M,
                                                A<:Array{T,3},
                                                O<:Poisson{T,AX,M,XYZ}
                                                }
    (;diag,axes,msk) = ope
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        r[i,j,k] = M(
            b[i,j,k]
            -ddi3(x,i,j,k)-ddj3(x,i,j,k)-ddk3(x,i,j,k)
                # -ddi(x,i,j,k,-1)-ddi(x,i,j,k,+1)
                # -ddj(x,i,j,k,-1)-ddj(x,i,j,k,+1)
                # -ddk(x,i,j,k,-1)-ddk(x,i,j,k,+1)
            +diag[i,j,k]*x[i,j,k],
            msk, i,j,k)
    end
end

function residual(r::A,x::A,b::A,ope::O) where {T,AX,M,
                                                A<:Array{T,3},
                                                O<:Poisson{T,AX,M,XY}
                                                }
    (;diag,axes,msk) = ope
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        r[i,j,k] = M(
            b[i,j,k]-ddi3(x,i,j,k)-ddj3(x,i,j,k)+diag[i,j,k]*x[i,j,k],
            msk, i,j,k)
    end
end


function jacobi(y::A,x::A,b::A,omega::T,ope::O) where {T,AX,M,
                                                       A<:Array{T,3},
                                                       O<:Poisson{T,AX,M,XYZ}}
    (;idiag,axes,msk)=ope
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        y[i,j,k] = M(
            (T(1)-omega)*x[i,j,k]-omega*idiag[i,j,k]*(
                b[i,j,k]-ddi3(x,i,j,k)-ddj3(x,i,j,k)-ddk3(x,i,j,k)
                # -ddi(x,i,j,k,-1)-ddi(x,i,j,k,+1)
                # -ddj(x,i,j,k,-1)-ddj(x,i,j,k,+1)
                # -ddk(x,i,j,k,-1)-ddk(x,i,j,k,+1)
                ),msk,i,j,k)
    end

end

function jacobi(y::A,x::A,b::A,omega::T,ope::O) where {T,AX,M,
                                                       A<:Array{T,3},
                                                       O<:Poisson{T,AX,M,XY}}
    (;idiag,axes,msk)=ope
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        y[i,j,k] = M(
            (T(1)-omega)*x[i,j,k]-omega*idiag[i,j,k]*(
                b[i,j,k]-ddi3(x,i,j,k)-ddj3(x,i,j,k)
                ),msk,i,j,k)
    end

end

function linerelaxation(x::A,b::A,rhs::V,d::V,ud::V,ope::O) where {T,AX,M,R,
                                                                   A<:Array{T,3},
                                                                   V<:Array{T,1},
                                                                   O<:Poisson{T,AX,M,R}}
    (;diag,axes,msk)=ope
    for j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        for k in CENTERS(axes.ax3)
            rhs[k] = M(b[i,j,k]-ddi3(x,i,j,k)-ddj3(x,i,j,k),
                       msk,i,j,k)
            d[k] = -diag[i,j,k]
            ud[k] = 1
        end
        tridiagsolve(view(x,i,j,:),rhs,d,ud)
    end

end

