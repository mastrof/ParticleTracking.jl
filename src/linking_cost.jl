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
    A = Blob(
        posA, CartesianIndex(ntuple(zero, N)),
        S(0 for _ in 1:N), zero(T), zero(T), zero(T), []
    )
    B = Blob(
        posB, CartesianIndex(ntuple(zero, N)),
        S(0 for _ in 1:N),
        zero(T), zero(T), zero(T), []
    )
    return cost(A, B; kwargs...)
end

"""
    QuadraticCost(A, B; g0=0, g2=0)
Evaluate the cost of linking two blobs `A` and `B` using the Euclidean
(𝐿²) norm.

The cost includes the spatial distance of the two blobs, and the distance
in their zeroth and second intensity moments.

**Keywords**
- `g0::Real = 0`: weight factor for the zeroth moment
- `g2::Real = 0`: weight factor for the second moment
"""
function QuadraticCost(A::AbstractBlob, B::AbstractBlob;
    g0::Real = 0.0, g2::Real = 0.0
)
    dx = location(B) .- location(A)
    dm0 = g0*(zeroth_moment(B) - zeroth_moment(A))
    dm2 = g2*(second_moment(B) - second_moment(A))
    return sqrt(sum(y*y for y in dx) + dm0*dm0 + dm2*dm2)
end

"""
    PCost(A, B; p=1, g0=0, g2=0)
Evaluate the cost of linking two blobs `A` and `B` using the 𝐿ᵖ norm.

The cost includes the spatial distance of the two blobs, and the distance
in their zeroth and second intensity moments.

**Keywords**
- `p::Real = 1`: the norm exponent; `p=1` corresponds to the Manhattan metric,
  `p=Inf` to the Chebyshev metric, and `p=2` to the Euclidean metric
  (in which case `QuadraticCost` provides slightly better performance)
- `g0::Real = 0`: weight factor for the zeroth moment
- `g2::Real = 0`: weight factor for the second moment
"""
function PCost(A::AbstractBlob, B::AbstractBlob;
    p::Real = 1, g0::Real = 0.0, g2::Real = 0.0
)
    dx = abs.(location(B) .- location(A))
    dm0 = g0*abs(zeroth_moment(B) - zeroth_moment(A))
    dm2 = g2*abs(second_moment(B) - second_moment(A))
    norm((dx..., dm0, dm2), p)
end
