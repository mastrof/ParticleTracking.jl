export connected_components, collect_groups

"""
    connected_components(stack::Stack; bkg = zero(eltype(first(stack))), dims = coords_spatial(first(stack)))
    connected_components(img::AbstractArray; bkg = zero(eltype(img)), dims = coords_spatial(img))
    connected_components(stack::Stack, connectivity; bkg = zero(eltype(first(stack))))
    connected_components(img::AbstractArray, connectivity; bkg = zero(eltype(img)))
Find the connected components in a stack `stack` or a single image `img`.
Components are defined as connected voxels that all have the same value,
distinct from the background component `bkg`.

Connectivity can be specified either through a list `dims` indicating which dimensions
are used to determine connectivity, or through a `connectivity` matrix.
See `label_components` from Images.jl for more details.
"""
function connected_components(stack::Stack;
    bkg = zero(eltype(first(stack))),
    dims = coords_spatial(first(stack))
)
    map(img -> connected_components(img; bkg, dims), stack)
end
function connected_components(img::AbstractArray;
    bkg = zero(eltype(img)), dims = coords_spatial(img)
)
    labels = label_components(img; bkg, dims)
    collect_groups(labels)
end

function connected_components(stack::Stack,
    connectivity; bkg = zero(eltype(first(stack)))
)
    map(img -> connected_components(img, connectivity; bkg), stack)
end
function connected_components(img::AbstractArray,
    connectivity; bkg = zero(eltype(img))
)
    labels = label_components(mg, connectivity; bkg)
    collect_groups(labels)
end


"""
    collect_groups(labels::AbstractArray{Int,D})
Group connected components labels. Works for `D=2` or `D=3`.

Each connected component is collected into a `Vector{CartesianIndex}`,
so the entire set of groups is returned as `Vector{Vector{CartesianIndex}}`.

`labels` can be obtained through the `label_components` function from Images.jl.
"""
function collect_groups(labels::AbstractMatrix{Int})
    groups = [CI{2}[] for _ in 1:maximum(labels)]
    for j in axes(labels,2), i in axes(labels,1)
        l = labels[i,j]
        if l ≠ 0
            push!(groups[l], CI(i,j))
        end
    end
    return groups
end
function collect_groups(labels::AbstractArray{Int,3})
    groups = [CI{3}[] for _ in 1:maximum(labels)]
    for k in axes(labels,3), j in axes(labels,2), i in axes(labels,1)
        l = labels[i,j,k]
        if l ≠ 0
            push!(groups[l], CI(i,j,k))
        end
    end
    return groups
end

