using BenchmarkTools

#include("axes.jl")
#include("boundary_conditions.jl")

# ICENTERS(a::Axis{HALOED,T}) where {T} = a.nhalo+2:a.nhalo+a.n-1
# ICENTERS(a::Axis{CLOSED,T}) where {T} = 2:a.n-1

coarse_to_fine(i, ax::Axis{HALOED,I}) where I = (i-ax.nhalo-1)*2+ax.nhalo+1
coarse_to_fine(i, ax::Axis{CLOSED,I}) where I = (i-1)*2+1

fine_to_coarse(i, ax::Axis{HALOED,I}) where I = div(i-ax.nhalo-1,2)+ax.nhalo+1
fine_to_coarse(i, ax::Axis{CLOSED,I}) where I = div(i-1,2)+1

# abstract type BoundaryCondition end
# struct Dirichlet<:BoundaryCondition end
# struct Neumann<:BoundaryCondition end

# abstract type Restriction end
# abstract type XYZ<:Restriction end
# abstract type XY<:Restriction end

abstract type RP{D,B} end


# function Restriction(D,B,axes)
#     coef = zeros(size(axes,CCC))
#     Restriction{D,B,eltype(coef),typeof(axes)}(coef,axes)
# end

# function Prolongation(D,B,axes)
#     coef = zeros(size(axes,CCC))
#     Prolongation{D,B,eltype(coef),typeof(axes)}(coef,axes)
# end

# (res::Restriction{D,B,T,AX})(xc,xf) where {D,B,T,AX} = restriction(xc,xf,res)


function restriction(xc::A,xf::A,coef::A,
                     axes::AX,::Type{R}) where {T,AX<:AXES,
                                                A<:Array{T,3},
                                                B<:BC,
                                                R<:RP{XYZ,B}}
    (;ax1,ax2,ax3) = axes

    @assert size(xc) == size(coef)

    for kk in CENTERS(ax3)
        k = coarse_to_fine(kk,ax3)
        for jj in CENTERS(ax2), ii in CENTERS(ax1)
            i = coarse_to_fine(ii,ax1)
            j = coarse_to_fine(jj,ax2)

            xc[ii,jj,kk] = coef[ii,jj,kk]*(
                xf[i,j  ,k  ] + xf[i+1,j  ,k  ]+
                    xf[i,j+1,k  ] + xf[i+1,j+1,k  ]+
                    xf[i,j  ,k+1] + xf[i+1,j  ,k+1]+
                    xf[i,j+1,k+1] + xf[i+1,j+1,k+1])

        end
    end
end

# function restriction(xc::A,xf::A,coef::A,
#                      axes::AX,::Type{R}) where {T,AX<:AXES,
#                                                 A<:Array{T,3},
#                                                 B<:BC,
#                                                 R<:RP{XY,B}}
#     (;ax1,ax2,ax3) = axes
#     for k in CENTERS(ax3)
#         for jj in CENTERS(ax2), ii in CENTERS(ax1)
#             i = coarse_to_fine(ii,ax1)
#             j = coarse_to_fine(jj,ax2)

#             xc[ii,jj,k] = coef[ii,jj,k]*(
#                 xf[i,j  ,k  ] + xf[i+1,j  ,k  ]+
#                     xf[i,j+1,k  ] + xf[i+1,j+1,k  ])

#         end
#     end
# end

function restriction(xc::A,xf::A,coef::A,
                     axes::AX,::Type{R}) where {T,AX<:AXES,
                                                A<:Array{T,3},
                                                B<:BC,
                                                R<:RP{XY,B}}
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(ax3)
        for jj in CENTERS(ax2), ii in CENTERS(ax1)
            i = coarse_to_fine(ii,ax1)
            j = coarse_to_fine(jj,ax2)

            xc[ii,jj,k] = coef[ii,jj,k]*Rsquare(xf,i,j,k,ax1,ax2)

        end
    end
end

@inline Rsquare(x,i,j,k,ax1,ax2) = (81(x[i,j,k]+x[i+1,j,k]+x[i,j+1,k]+x[i+1,j+1,k])
                            -9(ddi(x,i,j,k,-1,ax1)+ddi(x,i,j,k,2,ax1)
                               +ddj(x,i,j,k,-1)+ddj(x,i,j,k,2)
                               +ddij(x,i,j,k,-1,1)+ddij(x,i,j,k,2,1)
                               +ddij(x,i,j,k,1,-1)+ddij(x,i,j,k,1,2))
                            +(ddij(x,i,j,k,-1,-1)+ddij(x,i,j,k,-1,2)
                              +ddij(x,i,j,k,2,-1)+ddij(x,i,j,k,2,2))
                            )/64

#@inline square(x,i,j,k,di,dj) = 9x[i,j,k]+3x[i+di,j,k]+3x[i,j+dj,k]+x[i+di,j+dj,k]

@inline square(x,i,j,k,di,dj,ax1,ax2) = 9x[i,j,k]+3*ddi(x,i,j,k,di,ax1)+3*ddj(x,i,j,k,dj)+ddij(x,i,j,k,di,dj)

function prolongation(xf::A,xc::A,coef::A,
                     axes::AX,::Type{R}) where {T,AX<:AXES,
                                                A<:Array{T,3},
                                                B<:BC,
                                                R<:RP{XYZ,B}}

    @assert size(xf) == size(coef)

    (;ax1,ax2,ax3) = axes

    for kk in CENTERS(ax3)
        k = coarse_to_fine(kk,ax3)
        for jj in CENTERS(ax2), ii in CENTERS(ax1)
            i = coarse_to_fine(ii,ax1)
            j = coarse_to_fine(jj,ax2)

            ymm = square(xc,ii,jj,kk,-1,-1,ax1,ax2)
            ypm = square(xc,ii,jj,kk, 1,-1,ax1,ax2)
            ymp = square(xc,ii,jj,kk,-1, 1,ax1,ax2)
            ypp = square(xc,ii,jj,kk, 1, 1,ax1,ax2)

            if kk > 1
                xmm = square(xc,ii,jj,kk-1,-1,-1,ax1,ax2)
                xpm = square(xc,ii,jj,kk-1, 1,-1,ax1,ax2)
                xmp = square(xc,ii,jj,kk-1,-1, 1,ax1,ax2)
                xpp = square(xc,ii,jj,kk-1, 1, 1,ax1,ax2)
            else
                if B == DIRICHLET
                    xmm = -ymm
                    xmp = -ymp
                    xpm = -ypm
                    xpp = -ypp
                else
                    xmm = ymm
                    xmp = ymp
                    xpm = ypm
                    xpp = ypp
                end
            end
            if kk < ax3.n
                zmm = square(xc,ii,jj,kk+1,-1,-1,ax1,ax2)
                zpm = square(xc,ii,jj,kk+1, 1,-1,ax1,ax2)
                zmp = square(xc,ii,jj,kk+1,-1, 1,ax1,ax2)
                zpp = square(xc,ii,jj,kk+1, 1, 1,ax1,ax2)
            else
                if B == DIRICHLET
                    zmm = -ymm
                    zmp = -ymp
                    zpm = -ypm
                    zpp = -ypp
                else
                    zmm = ymm
                    zmp = ymp
                    zpm = ypm
                    zpp = ypp
                end
            end

            xf[i,j,k]     += coef[i,j,k]    *(3ymm+xmm)
            xf[i+1,j,k]   += coef[i+1,j,k]  *(3ypm+xpm)
            xf[i,j+1,k]   += coef[i,j+1,k]  *(3ymp+xmp)
            xf[i+1,j+1,k] += coef[i+1,j+1,k]*(3ypp+xpp)

            xf[i,j,k+1]     += coef[i,j,k+1]    *(3ymm+zmm)
            xf[i+1,j,k+1]   += coef[i+1,j,k+1]  *(3ypm+zpm)
            xf[i,j+1,k+1]   += coef[i,j+1,k+1]  *(3ymp+zmp)
            xf[i+1,j+1,k+1] += coef[i+1,j+1,k+1]*(3ypp+zpp)
        end

    end
end

function prolongation(xf::A,xc::A,coef::A,
                     axes::AX,::Type{R}) where {T,AX<:AXES,
                                                A<:Array{T,3},
                                                B<:BC,
                                                R<:RP{XY,B}}

    @assert size(xf) == size(coef)

    (;ax1,ax2,ax3) = axes

    for k in CENTERS(ax3)
        for jj in CENTERS(ax2), ii in CENTERS(ax1)
            i = coarse_to_fine(ii,ax1)
            j = coarse_to_fine(jj,ax2)

            ymm = square(xc,ii,jj,k,-1,-1,ax1,ax2)
            ypm = square(xc,ii,jj,k, 1,-1,ax1,ax2)
            ymp = square(xc,ii,jj,k,-1, 1,ax1,ax2)
            ypp = square(xc,ii,jj,k, 1, 1,ax1,ax2)

            xf[i,j,k]     += coef[i,j,k]    *ymm
            xf[i+1,j,k]   += coef[i+1,j,k]  *ypm
            xf[i,j+1,k]   += coef[i,j+1,k]  *ymp
            xf[i+1,j+1,k] += coef[i+1,j+1,k]*ypp

        end
    end
end

# function prolongation(xf::A,xc::A,coef::A,
#                      axes::AX,::Type{R}) where {T,AX<:AXES,
#                                                 A<:Array{T,3},
#                                                 B<:BC,
#                                                 R<:RP{XY,B}}

#     @assert size(xf) == size(coef)

#     (;ax1,ax2,ax3) = axes

#     for k in CENTERS(ax3)
#         for ii in CENTERS(ax1), jj in CENTERS(ax2)
#             i = coarse_to_fine(ii,ax1)
#             j = coarse_to_fine(jj,ax2)


#             xf[i,j,k]     += coef[i,j,k]    *xc[ii,jj,k]
#             xf[i+1,j,k]   += coef[i+1,j,k]  *xc[ii,jj,k]
#             xf[i,j+1,k]   += coef[i,j+1,k]  *xc[ii,jj,k]
#             xf[i+1,j+1,k] += coef[i+1,j+1,k]*xc[ii,jj,k]

#         end
#     end
# end



function test_RP()

    nx,ny,nz = 100,100,6
    nhalo = 1

    location = (CENTERS,CENTERS,CENTERS)



    axi = Axis(HALOED,nx,nhalo)
    axj = Axis(HALOED,ny,nhalo)
    axk = Axis(CLOSED,nz,nhalo)

    fine_axes = AXES(axi,axj,axk)

    xf,r,b=[zeros(size(fine_axes,location)) for _ in 1:3]

    axi = Axis(HALOED,div(nx,2),nhalo)
    axj = Axis(HALOED,div(ny,2),nhalo)
    axk = Axis(CLOSED,div(nz,2),nhalo)

    coarse_axes = AXES(axi,axj,axk)

    coarse_shape = size(coarse_axes,location)
    xc = zeros(coarse_shape)
    fcoef = xf*0
    ccoef = xc*0

    xc[2,2,1] = 1
    @. xf = 0
    @. fcoef = 1
    @. ccoef = 1

    rp = RP{XYZ,DIRICHLET}
    restriction(xc,xf,ccoef, coarse_axes, rp)
    prolongation(xf,xc,fcoef, coarse_axes, rp)
    #restriction(xc,xf,ccoef,coarse_axes,RESTRICTION{XYZ,DIRICHLET})
    #prolongation(xf,xc,fcoef,coarse_axes,PROLONGATION{XYZ,DIRICHLET})
    rp
end


function test_RP2D()

    nx,ny,nz = 100,100,1
    nhalo = 1

    location = (CENTERS,CENTERS,CENTERS)



    axi = Axis(HALOED,nx,nhalo)
    axj = Axis(HALOED,ny,nhalo)
    axk = Axis(CLOSED,nz,nhalo)

    fine_axes = AXES(axi,axj,axk)

    xf,r,b=[zeros(size(fine_axes,location)) for _ in 1:3]

    axi = Axis(HALOED,div(nx,2),nhalo)
    axj = Axis(HALOED,div(ny,2),nhalo)
    axk = Axis(CLOSED,nz,nhalo)

    coarse_axes = AXES(axi,axj,axk)

    coarse_shape = size(coarse_axes,location)
    xc = zeros(coarse_shape)
    fcoef = xf*0
    ccoef = xc*0

    xc[2,2,1] = 1
    @. xf = 0
    @. fcoef = 1
    @. ccoef = 1

    rp = RP{XY,DIRICHLET}
    restriction(xc,xf,ccoef, coarse_axes, rp)
    prolongation(xf,xc,fcoef, coarse_axes, rp)
    #restriction(xc,xf,ccoef,coarse_axes,RESTRICTION{XYZ,DIRICHLET})
    #prolongation(xf,xc,fcoef,coarse_axes,PROLONGATION{XYZ,DIRICHLET})
    rp
end
