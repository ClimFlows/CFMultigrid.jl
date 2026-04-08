struct PoissonNonUniform{T,AX,M,D} <: OPERATOR
    axes:: AX
    cxx::T
    cyy::T
    czz::T
    diag:: Array{T,3}
    idiag:: Array{T,3}
    msk::Array{Int8,3}
end

function PoissonNonUniform(axes::AX,::Type{B},::Type{D};msk=nothing,cxx=1.0,cyy=1.0,czz=1.0) where {AX<:AXES, D<:DIRECTIONS, B<:BC}
    T = Float64
    shape = size(axes,(CENTERS,CENTERS,CENTERS))
    diag, idiag = [zeros(T,shape) for _ in 1:2]
    m, status = check_msk(msk, shape)
    ope = PoissonNonUniform{T,AX,status,D}(axes,T(cxx),T(cyy),T(czz),diag,idiag,m)
    set_poisson_diag(ope,B)
    ope
end

function set_poisson_diag(ope::PoissonNonUniform{T,AX,M,D}, ::Type{B}) where {T,AX,M,D<:DIRECTIONS,B<:BC}
    (;diag,idiag,msk,cxx,cyy,czz) = ope
    shape = size(diag)
    x,r,b = [zeros(T,shape) for _ in 1:3]
    @. x = 1
    coef = diagcoef(shape)
    if coef==6
        @. diag =  2*(cxx+cyy+czz)
    elseif coef==4
        @assert shape[3] == 1
        @. diag =  2*(cxx+cyy)
    else
        @assert false
    end
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



function residual(r::A,x::A,b::A,ope::O) where {T,AX,M,D,
                                                A<:Array{T,3},
                                                O<:PoissonNonUniform{T,AX,M,D}
                                                }
    (;diag,axes,msk,cxx,cyy,czz) = ope
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        r[i,j,k] =M(
        b[i,j,k]-(
            +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
            +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy
            +(ddk(x,i,j,k,-1,ax3)+ddk(x,i,j,k,+1,ax3))*czz
            -diag[i,j,k]*x[i,j,k]), msk, i,j,k)
    end
end

function jacobi(y::A,x::A,b::A,omega::T,ope::O) where {T,AX,M,D,
                                                       A<:Array{T,3},
                                                       O<:PoissonNonUniform{T,AX,M,D}}
    (;idiag,axes,msk,cxx,cyy,czz)=ope
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        y[i,j,k] = M(
            (T(1)-omega)*x[i,j,k]-omega*idiag[i,j,k]*(
                b[i,j,k]-(
                    +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
                    +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy
                    +(ddk(x,i,j,k,-1,ax3)+ddk(x,i,j,k,+1,ax3))*czz
                )),msk,i,j,k)
    end

end

function linerelaxation(x::A,b::A,rhs::V,d::V,ud::V,ope::O) where {T,AX,M,D,
                                                                   A<:Array{T,3},
                                                                   V<:Array{T,1},
                                                                   O<:PoissonNonUniform{T,AX,M,D}}
    (;diag,axes,msk,cxx,cyy,czz)=ope
    (;ax1,ax2,ax3) = axes
    for j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        for k in CENTERS(axes.ax3)
            rhs[k] = M(b[i,j,k]-(
                +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
                +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy),msk,i,j,k)
            d[k] = -diag[i,j,k]
            ud[k] = czz
        end
        tridiagsolve(view(x,i,j,:),rhs,d,ud)
    end

end
