"""

 Helmholtz operator

The `diag` term is `diag = C+D` with
  - D = 1/2(dx^2+dy^2+dz^2), the Poisson part
  - C, the identity part

`idiag` is `idiag = 1/diag`, like for Poisson

During the coarsening, we keep D constant (it should be divided by
4). To compensate, C must multiplied by 4.

During the coarsening, `restriction` averages x then multiplies by
4. To set the diagonal on the coarsened grid the technique is simple:
apply the restriction on the coefficient. This is one in `set_ope_coef`

"""
struct Helmholtz{T,AX,M,D} <: OPERATOR
    axes:: AX
    cxx::T
    cyy::T
    czz::T
    coef:: Array{T,3}
    diag:: Array{T,3}
    idiag:: Array{T,3}
    msk::Array{Int8,3}
end

function Helmholtz(axes::AX,::Type{B},::Type{D};msk=nothing,cxx=1.0,cyy=1.0,czz=1.0,coef=1.0) where {AX<:AXES, D<:DIRECTIONS, B<:BC}
    T = Float64
    shape = size(axes,(CENTERS,CENTERS,CENTERS))
    diag, idiag,identitycoef = [zeros(T,shape) for _ in 1:3]
    m, status = check_msk(msk, shape)
    @. identitycoef = coef
    ope = Helmholtz{T,AX,status,D}(axes,T(cxx),T(cyy),T(czz),identitycoef,diag,idiag,m)
    set_helmholtz_diag(ope,B)
    ope
end

get_RP(grid::G) where {RP,AX,T,O,G<:Grid{RP,AX,T,O}} = RP

function set_ope_coef(fine,coarse,::Type{H},::Type{B}) where {T,AX,M,D<:DIRECTIONS,B<:BC,H<:Helmholtz{T,AX,M,D}}
    RP = get_RP(coarse)
    restriction(coarse.ope.coef,fine.ope.coef,coarse.Rcoef,coarse.axes,RP)
    set_helmholtz_diag(coarse.ope, B)
end


function set_helmholtz_diag(ope::Helmholtz{T,AX,M,D}, ::Type{B}) where {T,AX,M,D<:DIRECTIONS,B<:BC}
    (;diag,idiag,msk,cxx,cyy,czz,coef) = ope
    shape = size(diag)
    x,r,b = [zeros(T,shape) for _ in 1:3]
    @. x = 1
    ncoef = diagcoef(shape)
    if ncoef==6
        @. diag =  2*(cxx+cyy+czz)
    elseif ncoef==4
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
    # we add the identity term here
    @. diag = diag + coef


    if M==ON
        @. idiag[msk==1]  =  1/diag[msk==1]
    else
        @. idiag  =  1/diag
    end
end



function residual(r::A,x::A,b::A,ope::O) where {T,AX,M,
                                                A<:Array{T,3},
                                                O<:Helmholtz{T,AX,M,XYZ}
                                                }
    (;diag,axes,msk,cxx,cyy,czz) = ope
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        r[i,j,k] =M(
        b[i,j,k]-(
            # +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
            # +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy
            # +(ddk(x,i,j,k,-1,ax3)+ddk(x,i,j,k,+1,ax3))*czz
            ddi3(x,i,j,k)*cxx+ddj3(x,i,j,k)*cyy+ddk3(x,i,j,k)*czz
            -diag[i,j,k]*x[i,j,k]), msk, i,j,k)
    end
end

function residual(r::A,x::A,b::A,ope::O) where {T,AX,M,
                                                A<:Array{T,3},
                                                O<:Helmholtz{T,AX,M,XY}
                                                }
    (;diag,axes,msk,cxx,cyy,czz) = ope
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        r[i,j,k] =M(
        b[i,j,k]-(
            # +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
            # +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy
            ddi3(x,i,j,k)*cxx+ddj3(x,i,j,k)*cyy
            -diag[i,j,k]*x[i,j,k]), msk, i,j,k)
    end
end

function jacobi(y::A,x::A,b::A,omega::T,ope::O) where {T,AX,M,
                                                       A<:Array{T,3},
                                                       O<:Helmholtz{T,AX,M,XYZ}}
    (;idiag,axes,msk,cxx,cyy,czz)=ope
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        y[i,j,k] = M(
            (T(1)-omega)*x[i,j,k]-omega*idiag[i,j,k]*(
                b[i,j,k]-(
                    # +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
                    # +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy
                    # +(ddk(x,i,j,k,-1,ax3)+ddk(x,i,j,k,+1,ax3))*czz
                    ddi3(x,i,j,k)*cxx+ddj3(x,i,j,k)*cyy+ddk3(x,i,j,k)*czz
                )),msk,i,j,k)
    end

end

function jacobi(y::A,x::A,b::A,omega::T,ope::O) where {T,AX,M,
                                                       A<:Array{T,3},
                                                       O<:Helmholtz{T,AX,M,XY}}
    (;idiag,axes,msk,cxx,cyy,czz)=ope
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(axes.ax3), j in CENTERS(axes.ax2), i in CENTERS(axes.ax1)
        y[i,j,k] = M(
            (T(1)-omega)*x[i,j,k]-omega*idiag[i,j,k]*(
                b[i,j,k]-(
                    # +(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,+1,ax1))*cxx
                    # +(ddj(x,i,j,k,-1,ax2)+ddj(x,i,j,k,+1,ax2))*cyy
                    ddi3(x,i,j,k)*cxx+ddj3(x,i,j,k)*cyy
                )),msk,i,j,k)
    end

end

function linerelaxation(x::A,b::A,rhs::V,d::V,ud::V,ope::O) where {T,AX,M,D,
                                                                   A<:Array{T,3},
                                                                   V<:Array{T,1},
                                                                   O<:Helmholtz{T,AX,M,D}}
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
