import Base.@kwdef
using ManagedLoops

ie(grid) = return 1+grid.nh:grid.nx+grid.nh+1
je(grid) = return 1+grid.nh:grid.ny+grid.nh+1
ke(grid) = return 1+grid.nh:grid.nz+grid.nh+1

iv(grid) = return 1+grid.nh+1:grid.nx+grid.nh
jv(grid) = return 1+grid.nh+1:grid.ny+grid.nh
kv(grid) = return 1+grid.nh+1:grid.nz+grid.nh

ic(grid) = return 1+grid.nh:grid.nx+grid.nh
jc(grid) = return 1+grid.nh:grid.ny+grid.nh
kc(grid) = return (grid.case == :twod ? (1:grid.nz) : (1+grid.nh:grid.nz+grid.nh))

icp(grid) = return grid.nh:grid.nx+grid.nh+1
jcp(grid) = return grid.nh:grid.ny+grid.nh+1
kcp(grid) = return grid.nh:grid.nz+grid.nh+1


ia(grid) = return 1:grid.tx
ja(grid) = return 1:grid.ty
ka(grid) = return 1:grid.tz

c2fidx(i, nh) = return (i-nh-1)*2+nh+1
f2cidx(i, nh) = return div(i-nh-1,2)+nh+1

# c2fidx(i, nh) = return (i-nh)*2+nh
# f2cidx(i, nh) = return div(i-nh,2)+nh

A = Array{Float64, 3}
M = Array{UInt8, 3}

"""
    param = Param()
returns the default parameters set for the multigrid. `param` is mutable. All the entries can be modified.

    print(param)
to see all the parameters. `param` is the sole input argument to create a multigrid object (see ?get_gmg)
"""
@kwdef mutable struct Param
    nx = 8*8
    ny = 8*8
    nz = 8*8
    nh = 3
    npre = 1
    npost = 1
    ndeepest = 20
    maxite = 15
    tol = 1e-10
    omega = 0.9
#    maxlevs = 99
    location = :centers
    case = :threed
    operator = :laplacian
    cycle = :v
    verbose = true
    topology = :closed
end

Base.show(io::IO, p:: Param) = begin
    println("geometric multigrid parameters:")
    for name in propertynames(p)
        println(io, "  - ", name, ": ", getfield(p, name))
    end
end

struct Grid
    nx :: Int
    ny :: Int
    nz :: Int
    nh :: Int
    tx :: Int
    ty :: Int
    tz :: Int
    case:: Symbol
    function Grid(nx,ny,nz,nh,case)
        @assert case in [:twod, :threed] "grid must be :twod or :threed"
        tx = nx+2*nh
        ty = ny+2*nh
        case == :twod ? tz = nz : tz=nz+2*nh
        new(nx,ny,nz,nh,tx,ty,tz,case)
    end
end

struct Data
    x:: Array{Float64, 3}
    r:: Array{Float64, 3}
    b:: Array{Float64, 3}
    y:: Array{Float64, 3}
    function Data(grid:: Grid)
        x = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        r = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        b = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        y = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        new(x,r,b,y)
    end
end

struct Level
    grid:: Grid
end

struct Operators
    diag:: Array{Float64, 3}
    idiag:: Array{Float64, 3}
    Rcoef:: Array{Float64, 3}
    Pcoef:: Array{Float64, 3}
    msk:: Array{UInt8, 3}
    function Operators(grid:: Grid)
        diag = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        idiag = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        Rcoef = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        Pcoef = zeros(Float64, (grid.tx, grid.ty, grid.tz))
        msk = zeros(UInt8, (grid.tx, grid.ty, grid.tz))
        new(diag,idiag,Rcoef,Pcoef,msk)
    end
end


struct Gmg
    mgr:: Vector{ManagedLoops.HostManager}
    nlevels:: Integer
    param:: Param
    levels:: Vector{Level}
    data:: Vector{Data}
    oper:: Vector{Operators}
    smoother:: Vector{Function}
    residual:: Vector{Function}
    restriction:: Vector{Function}
    prolongation:: Vector{Function}
end
