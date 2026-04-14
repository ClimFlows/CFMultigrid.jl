#include("axes.jl")
#include("boundary_conditions.jl")

ICENTERS(a::Axis{HALOED,T}) where {T} = a.nhalo+2:a.nhalo+a.n-1
ICENTERS(a::Axis{CLOSED,T}) where {T} = 2:a.n-1


#using Plots

struct Grid{RP,AX,T,O}
    axes:: AX
    x:: Array{T,3}
    b:: Array{T,3}
    r:: Array{T,3}
    y:: Array{T,3}
    Rcoef:: Array{T,3}
    Pcoef:: Array{T,3}
    ope:: O
end

Operator(grid::G) where {RP,AX,T,O,G<:Grid{RP,AX,T,O}} = O

function Grid(axes,OPE,bc,switch;kwargs...)
    x,b,r,y,Rcoef,Pcoef=[zeros(size(axes,CCC)) for _ in 1:6]
    ope = OPE(axes,bc,switch;kwargs...)
    rp = RP{switch,bc}
    AX=typeof(axes)
    T=eltype(x)
    O=typeof(ope)
    Grid{rp,AX,T,O}(axes,x,b,r,y,Rcoef,Pcoef,ope)
end

# Base.getindex(x::Array{T,3},i::I,j::I,k::I) where {T,I<:Int64} = begin
#     if (0<i<=size(x,1)) & (0<j<=size(x,2)) & (0<k<=size(x,3))
#         return x[i,j,k]
#     else
#         return T(0)
#     end
# end





#include("lowlevel_operators.jl")
#include("RPs.jl")
#include("smoothers.jl")



function test_grid()
    nx,ny,nz = 30,30,1
    nhalo = 1


    axi = Axis(CLOSED,nx,nhalo)
    axj = Axis(CLOSED,ny,nhalo)
    axk = Axis(CLOSED,nz,nhalo)

    axes = AXES(axi,axj,axk)

    msk = zeros(Int8,size(axes,(CENTERS,CENTERS,CENTERS)))
    @. msk = 1
    @. msk[end-div(nx,2):end,1:10,:] = 0


    bc = DIRICHLET
    smoother = Jacobi(0.9)
    grid = Grid(axes,Poisson,bc,XYZ;msk=msk)
    (;b,x,r,ope) = grid


    b[1+div(nx,3),1+div(ny,3),1] = 1
    b[nx-div(nx,3),ny-div(ny,3),1] = -1
    @. x = 0

    smooth(1000, grid, smoother)

    residual(r,x,b,ope)

    print(maximum(@. abs(r)))
    levels = Vector(-.5:.05:.5)
    contourf(x[:,:,1];cmap=:RdBu)

end
