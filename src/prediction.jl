struct PredictiveCost{F, V<:Tuple}
    cost_func::F
    vel::Dict{Tuple{Int,Int}, V}
    weight::Float64
end
function (pc::PredictiveCost)(A, B, from_node::Tuple{Int,Int}, dt::Integer; kwargs...)
    # base spatial cost
    cost = pc.cost_func(A, B; kwargs...)
    # predict next position
    velocity = get(pc.vel, from_node, nothing)
    if isnothing(velocity) # cannot predict
        return cost
    end
    predicted_pos = location(A) .+ velocity .* dt
    prediction_cost = sum(abs2.(location(B) .- predicted_pos))
    # weighted cost
    (1-pc.weight)*cost + pc.weight*prediction_cost
end

function estimate_velocities!(
    velocities::AbstractDict,
    ancestor_map::AbstractDict,
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob{T,N}}},
    t::Integer,
    k::Integer
) where {T,N}
    k == 0 && return
    v = zeros(N)
    for idx in eachindex(blobs[t])
        # to evaluate v there must be a parent node
        current_node = get(ancestor_map, (t, idx), nothing)
        isnothing(current_node) && continue
        # preallocate k parents
        # but later calculation will stop at the last real parent
        # this improves performance by avoiding reallocations
        ancestors = fill(current_node, k)
        # find parent nodes in the last k steps
        num_parents = 1
        for step in 1:k-1
            parent_node = get(ancestor_map, current_node, nothing)
            if isnothing(parent_node)
                break
            end
            ancestors[step+1] = parent_node
            num_parents += 1
            current_node = parent_node
        end
        fill!(v, 0) # reset v
        for n in 1:num_parents # larger n => further in the past
            t2, idx2 = n == 1 ? (t, idx) : ancestors[n-1]
            t1, idx1 = ancestors[n]
            pos2 = location(blobs[t2][idx2])
            pos1 = location(blobs[t1][idx1])
            dt = t2 - t1
            @. v += (pos2 - pos1) / dt
        end
        v ./= num_parents
        velocities[(t, idx)] = Tuple(v)
    end
end
