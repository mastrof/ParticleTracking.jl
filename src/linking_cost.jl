export QuadraticCost, PCost, evaluate_costs!

"""
    evaluate_costs!(C, from, to, dt, cost, maxdist; kwargs...)
Fill the cost matrix `C` with the costs for linking two vectors of blobs
(`from` and `to`) `dt` frames apart.

The cost is evaluated according to the `cost` function and
`maxdist` sets the maximal distance (per frame) allowed for the link
between two blobs.

Optional kwargs can be passed to the `cost` function.
"""
function evaluate_costs!(C::OffsetMatrix, 
    from::AbstractVector{<:AbstractBlob},
    to::AbstractVector{<:AbstractBlob},
    dt::Integer,
    cost::Function,
    maxdist::Real;
    kwargs...
)
    for j in axes(C,2), i in axes(C,1)
        if i != 0 && j != 0
            C[i,j] = cost(from[i], to[j]; kwargs...)
            # costs larger than those for creation/annihilation
            # are set to Inf to improve performance
            C[i,j] > (dt*maxdist)^2 && (C[i,j] = Inf)
        elseif xor(i == 0, j == 0)
            C[i,j] = (dt*maxdist)^2
        end
    end
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
    dx = sqrt(sum(abs2.(xB .- xA)))
    dm0 = abs(zeroth_moment(B) - zeroth_moment(A))
    dm2 = abs(second_moment(B) - second_moment(A))
    return dx^p + g0*dm0^p + g2*dm2^p
end
