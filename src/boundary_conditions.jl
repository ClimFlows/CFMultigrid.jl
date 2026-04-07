include("axes.jl")

abstract type BC end
abstract type NONPERIO<:BC end
struct NEUMANN<:NONPERIO end
struct DIRICHLET<:NONPERIO end
struct PERIO<:BC end

abstract type BC_METHOD end

struct SINGLE{AX,BCS}<:BC_METHOD
    axes:: AX
    bcs:: BCS
end

function apply_bc(x, method::BM) where {BM<:SINGLE}
    (xbc,ybc,zbc) = method.bcs
    (;ax1,ax2,ax3) = method.axes

    if hashalo(ax1)
        @. x[1,:,:] = xbc(x[2,:,:], x[end-1,:,:])
        @. x[end,:,:] = xbc(x[end-1,:,:], x[2,:,:])
    end

    if hashalo(ax2)
        @. x[:,1,:] = ybc(x[:,2,:], x[:,end-1,:])
        @. x[:,end,:] = ybc(x[:,end-1,:], x[:,2,:])
    end

end

# function lateral_bc(x::A, axes::AX, ::Type{B}) where {T,
#                                                       AX<:AXES,
#                                                       A<:Array{T,3},
#                                                       B<:Boundary{NONPERIO,NONPERIO,NONPERIO}}
# end

hashalo(ax::AX) where {A,T,AX<:Axis{A,T}} = (A==HALOED)

# function fillhalo(x::A, axes::AX, ::Type{B}) where {T,
#                                                       AX<:AXES,
#                                                       A<:Array{T,3},
#                                                       B<:Boundary{PERIO,NONPERIO,NONPERIO}}
#     (;ax1) = axes
#     n = ax1.n
#     @assert hashalo(ax1)
#     if ax1.ntiles == 1
#         for k in CENTERS(ax3), j in CENTERS(ax2)
#             x[1,j,k] = x[1+n,j,k]
#             x[n+2,j,k] = x[2,j,k]
#         end
#     end
# end


(::Type{DIRICHLET})(xnext, xopposite) = -xnext
(::Type{NEUMANN})(xnext, xopposite) = xnext
(::Type{PERIO})(xnext, xopposite) = xopposite



function apply_lateral_bc(x::A,
                          axes::AX,
                          ::Type{B}) where {T,
                                            AX<:AXES,
                                            A<:Array{T,3},
                                            B<:BC}
    (;ax1,ax2,ax3) = axes
    for k in CENTERS(ax3)
        x[1,1,k] = B(x[2,2,k])
        x[end,end,k] = B(x[end-1,end-1,k])
        x[1,end,k] = B(x[2,end-1,k])
        x[end,1,k] = B(x[end-1,2,k])

        x[1,2:end-1,k] = B(x[2,2:end-1,k])
        x[2:end-1,1,k] = B(x[2:end-1,2,k])
        x[end,2:end-1,k] = B(x[end-1,2:end-1,k])
        x[2:end-1,end,k] = B(x[2:end-1,end-1,k])
    end
end

nx,ny,nz,nhalo = 4,6,1,1
location = (HALOED,HALOED,CLOSED)
axes = Axes(nx,ny,nz,nhalo,location...)
shape = size(axes,CCC)

x = zeros(shape)
for i in eachindex(x)
    x[i]=i
end

bcmeth = SINGLE(axes,(PERIO,PERIO,DIRICHLET))

apply_bc(x, bcmeth)
