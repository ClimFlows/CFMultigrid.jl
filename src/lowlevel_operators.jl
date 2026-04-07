abstract type OPERATOR end

@inline ddi(x::A,i,j,k,di,ax::Axis{HALOED,I}) where {I,T,A<:Array{T,3}} = x[i+di,j,k]
@inline ddi(x::A,i,j,k,di,ax::Axis{CLOSED,I}) where {I,T,A<:Array{T,3}} = (0<i+di<=ax.n) ? x[i+di,j,k] : T(0)

# @inline function ddi(x::A,i,j,k,di) where{T,A<:Array{T,3}}
#     if 0<i+di<=size(x,1)
#         return x[i+di,j,k]
#     else
#         return T(0)
#     end
# end

@inline function ddj(x::A,i,j,k,dj) where{T,A<:Array{T,3}}
    if 0<j+dj<=size(x,2)
        return x[i,j+dj,k]
    else
        return T(0)
    end
end

@inline function ddk(x::A,i,j,k,dk) where{T,A<:Array{T,3}}
    if 0<k+dk<=size(x,3)
        return x[i,j,k+dk]
    else
        return T(0)
    end
end

@inline function ddij(x::A,i,j,k,di,dj) where{T,A<:Array{T,3}}
    if (0<j+dj<=size(x,2)) && (0<i+di<size(x,1))
        return x[i+di,j+dj,k]
    else
        return T(0)
    end
end

@inline function ddi3(x::A,i,j,k) where{T,A<:Array{T,3}}
    if i==1
        return x[i+1,j,k]
    elseif i<size(x,1)
        return x[i-1,j,k]+x[i+1,j,k]
    else
        return x[i-1,j,k]
    end
end
@inline function ddj3(x::A,i,j,k) where{T,A<:Array{T,3}}
    if j==1
        return x[i,j+1,k]
    elseif j<size(x,2)
        return x[i,j-1,k]+x[i,j+1,k]
    else
        return x[i,j-1,k]
    end
end
@inline function ddk3(x::A,i,j,k) where{T,A<:Array{T,3}}
    if k==1
        return x[i,j,k+1]
    elseif k<size(x,3)
        return x[i,j,k-1]+x[i,j,k+1]
    else
        return x[i,j,k-1]
    end
end



# @inline ddi(x::A,i,j,k,di,ax::Axis{HALOED,I}) where {I,T,A<:Array{T,3}} = x[i+di,j,k]
# @inline ddi(x::A,i,j,k,di,ax::Axis{CLOSED,I}) where {I,T,A<:Array{T,3}} = (0<i+di<=ax.n) ? x[i+di,j,k] : T(0)
# @inline ddi(x::A,i,j,k,di,ax::Axis{EMPTY,I}) where {I,T,A<:Array{T,3}} = T(0)

# @inline ip(x,i,j,k,axes) = ddi(x,i,j,k,1,axes.ax1)


abstract type MASKED end
abstract type OFF<:MASKED end
abstract type ON<:MASKED end

@inline (::Type{ON})(x,msk,i,j,k) = msk[i,j,k]*x
@inline (::Type{OFF})(x,msk,i,j,k) = x


diagcoef(shape) = 2*sum(s>1 for s in shape)


function check_msk(msk,shape)
    if msk == nothing
        return zeros(Int8,1,1,1), OFF
    else
        @assert eltype(msk) == Int8
        @assert size(msk) == shape
        return msk, ON
    end
end


include("poisson.jl")
include("poisson_nonuniform.jl")

function set_ope_coef(fine,coarse) end



function tridiagsolve(x,rhs,d,ud)
    n = length(x)
    for k in 2:n
        w = ud[k]/d[k-1]
        d[k] -= w*ud[k-1]
        rhs[k] -= w*rhs[k-1]
    end
    x[n] = rhs[n]/d[n]
    for k in n-1:-1:1
        x[k] = (rhs[k]-ud[k]*x[k+1])/d[k]
    end
end

function test_tridiag()
    n = 10
    x, rhs,d,ud,y,b0,d0 = [zeros(n) for _ in 1:7]

    rhs[1] = 1
    #rhs = rand(n)
    @. d = 2
    @. ud = -1

    @. b0 = rhs
    @. d0 = d

    tridiagsolve(x,rhs,d,ud)
    #println(x)

    @. y = d0*x
    for i in 1:n
        if i>1
            y[i]+=ud[i]*x[i-1]
        end
        if i<n
            y[i]+=ud[i]*x[i+1]
        end
    end
    #println(y-b0)

    @assert all(isapprox.(y-b0,0,atol=1e-12))
end
