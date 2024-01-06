export QuadraticCost, PCost

"""
    evaluate_costs!(C, from, to, dt, cost, maxcost; kwargs...)
Fill the cost matrix `C` with the costs for linking two vectors of blobs
(`from` and `to`) `dt` frames apart.

Optional kwargs can be passed to the `cost` function.
"""
function evaluate_costs!(C::OffsetMatrix, 
    from::AbstractVector{T},
    to::AbstractVector{T},
    dt::Integer,
    cost::Function,
    maxcost::Real;
    kwargs...
) where {T<:AbstractBlob}
    for j in axes(C,2), i in axes(C,1)
        if i != 0 && j != 0
            C[i,j] = cost(from[i], to[j]; kwargs...)
            # costs larger than those for creation/annihilation
            # are set to Inf to improve performance
            if C[i,j] > maxcost
                C[i,j] = Inf
            end
        elseif xor(i == 0, j == 0)
            C[i,j] = maxcost
        end
    end
end

"""
    evaluate_maxcost(T, maxdist, dt, cost; kwargs...)
Evaluate the maximum allowed cost for a link if the maximum allowed
distance is `maxdist` pixels and the time difference `dt` frames.
"""
function evaluate_maxcost(::Type{<:AbstractBlob{T,S,N}},
    maxdist::Real, dt::Integer, cost::Function; kwargs...
) where {T,S,N}
    # define two dummy blobs dt*maxdist apart and evaluate their linking cost
    posA = ntuple(_ -> 0.0, N)
    posB = ntuple(i -> i==1 ? Float64(dt*maxdist) : 0.0, N)
    A = BlobRefined(
        posA, CartesianIndex(ntuple(zero, N)),
        S(0 for _ in 1:N), zero(T), zero(T), zero(T), []
    )
    B = BlobRefined(
        posB, CartesianIndex(ntuple(zero, N)),
        S(0 for _ in 1:N),
        zero(T), zero(T), zero(T), []
    )
    return cost(A, B; kwargs...)
end

function QuadraticCost(A::AbstractBlob, B::AbstractBlob;
    g0::Real = 1.0, g2::Real = 1.0
)
    xA = location(A)
    xB = location(B)
    dx2 = sum(abs2.(xB .- xA))
    dm02 = (zeroth_moment(B) - zeroth_moment(A))^2
    dm22 = (second_moment(B) - second_moment(A))^2
    return dx2 + g0*dm02 + g2*dm22
end

function PCost(A::AbstractBlob, B::AbstractBlob;
    p::Real = 1, g0::Real = 1.0, g2::Real = 1.0
)
    xA = location(A)
    xB = location(B)
    dxp = sum(@. abs(xB - xA)^p)
    dm0p = abs(zeroth_moment(B) - zeroth_moment(A))^p
    dm2p = abs(second_moment(B) - second_moment(A))^p
    return dxp + g0*dm0p + g2*dm2p
end
