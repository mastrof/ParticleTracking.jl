function initialize_links!(G::OffsetMatrix, C::OffsetMatrix,
    from::AbstractVector{<:AbstractBlob}, to::AbstractVector{<:AbstractBlob},
    t::Integer, r::Integer
)
    js_occupied = zero(axes(G, 2))
    # initial guesses for min cost links
    for i in axes(G,1)[1:end]
        # TODO: add check to avoid throwing errors when empty
        c, j = findmin(C[i,j] for j in axes(G,2)[1:end])
        if js_occupied[j] == 0 && !isinf(C[i,j])
            G[i,j] = 1
            if j != 0
                js_occupied[j] = 1
            end
        else
            G[i,0] = 1
        end
    end
    # link all unassigned blobs to dummies
    for j in axes(G,2)[1:end]
        if sum(G[:,j]) == 0
            G[0,j] = 1
        end
    end
end


function optimize_links!(G::AbstractMatrix, C::AbstractMatrix)
    Z = similar(C, Float64)
    fill!(Z, -Inf)
    while any(<(0), Z)
        optimize_step!(G, Z, C)
    end
end

function optimize_step!(
    G::AbstractMatrix, Z::AbstractMatrix, C::AbstractMatrix
)
    fill!(Z, 0.0) # Reset reduced costs
    # look for non-linked pairs with finite cost
    indices = findall(idx -> G[idx] == 0 && !isinf(C[idx]), CartesianIndices(G))
    minZ = 0.0
    minIJKL = [0, 0, 0, 0]
    for (I, J) in Tuple.(indices)
        if I != 0 && J != 0 # swap
            L = findfirst(G[I,j] == 1 for j in axes(G,2))
            K = findfirst(G[i,J] == 1 for i in axes(G,1))
            Z[I,J] = C[I,J] - C[I,L] - C[K,J] + C[K,L]
        elseif I == 0 && J > 0 # creation
            L = 0
            K = findfirst(G[i,J] == 1 && i > 0 for i in axes(G,1))
            Z[0,J] = C[0,J] - C[K,J] + C[K,0]
        elseif I > 0 && J == 0 # annihilation
            L = findfirst(G[I,j] == 1 && j > 0 for j in axes(G,2))
            K = 0
            Z[I,0] = C[I,0] - C[I,L] + C[0,L]
        end
        if Z[I, J] < minZ
            minZ = Z[I, J]
            minIJKL .= (I, J, K, L)
        end
    end
    # update links with minimal reduced cost if any
    if minZ < 0
        I, J, K, L = minIJKL
        G[I, J] = 1
        if I != 0
            G[I, L] = 0
        end
        if J != 0
            G[K, J] = 0
        end
        G[K, L] = 1
    end
    return nothing
end

function makelinks(
    linking_matrix::OffsetMatrix,
    from_indices::AbstractVector{<:Integer},
    to_indices::AbstractVector{<:Integer},
    t::Integer,
    r::Integer
)
    links = NamedTuple{(:from_t, :from_idx, :to_t, :to_idx), NTuple{4,Int}}[]
    for i in eachindex(from_indices)
        row = view(linking_matrix, i, :)
        j = findfirst(==(1), row)
        if !isnothing(j) && j > 0 # real link between two blobs
            link = (
                from_t=t,
                from_idx=from_indices[i],
                to_t=t+r,
                to_idx=to_indices[j],
            )
            push!(links, link)
        end
    end
    links
end
