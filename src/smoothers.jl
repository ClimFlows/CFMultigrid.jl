abstract type SMOOTHER end

struct Jacobi{T}<:SMOOTHER
    omega::T
end

struct Gauss_Seidel<:SMOOTHER end

struct LineRelaxation{T}<:SMOOTHER
    rhs::Array{T,1}
    d::Array{T,1}
    ud::Array{T,1}
end

function LinearRelaxation(n,T)
    rhs,d,ud = [zeros(T,n) for _ in 1:3]
    LineRelaxation{T}(rhs,d,ud)
end

function smooth(n, grid, smoother::Jacobi)
    (;x,y,b,ope) = grid
    (;omega)=smoother
    for _ in 1:n
        jacobi(y,x,b,omega,ope)
        jacobi(x,y,b,omega,ope)
    end
end

function smooth(n, grid, smoother::Gauss_Seidel)
    (;x,b,ope) = grid
    for _ in 1:n
        gauss_seidel(x,b,ope)
    end
end

function smooth(n, grid, smoother::LineRelaxation)
    (;x,b,ope) = grid
    (;rhs,d,ud) = smoother
    for _ in 1:n
        linerelaxation(x,b,rhs,d,ud,ope)
    end
end

#(smoother::S)(n,grid) where {S<:SMOOTHER} = smooth(n,grid,smoother)
