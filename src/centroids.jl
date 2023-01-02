export component_centroids

"""
    component_centroids(c::AbstractVector{<:AbstractVector{<:CC}})
Evaluate centroids of connected components `c`
(a vector of connected components for each frame in a stack)
in image coordinates.
"""
Images.component_centroids(c::AbstractVector{<:AbstractVector{<:CC}}) = component_centroids.(c)
"""
    component_centroids(c::AbstractVector{<:CC})
Evaluate centroids of connected components `c`
(a vector of connected components) in image coordinates.
"""
Images.component_centroids(c::AbstractVector{<:CC}) = component_centroids.(c)
"""
    component_centroids(c::CC)
Evaluate centroid of a connected component `c` in image coordinates.
"""
Images.component_centroids(c::CC) = Tuple(sum(c)) ./ length(c)


"""
    component_centroids(c::AbstractVector{<:AbstractVector{CC}}, coords)
Evaluate centroids of connected components `c`
(a vector of connected components for each frame in a stack)
in real-space coordinates `coords`.

`coords` should be an iterable where each element `coords[i]` is the
mesh associated to the `i`-th spatial dimension of the original image.
"""
function Images.component_centroids(c::AbstractVector{<:AbstractVector{<:CC}}, coords)
    map(p -> component_centroids(p, coords), c)
end
"""
    component_centroids(c::AbstractVector{CC}, coords)
Evaluate centroids of connected components `c`
(a vector of connected components) in real-space coordinates `coords`.

`coords` should be an iterable where each element `coords[i]` is the
mesh associated to the `i`-th spatial dimension of the original image.
"""
function Images.component_centroids(c::AbstractVector{<:CC}, coords)
    map(p -> component_centroids(p, coords), c)
end
"""
    component_centroids(c::CC, coords)
Evaluate centroid of a connected component `c` in real-space coordinates `coords`.

`coords` should be an iterable where each element `coords[i]` is the
mesh associated to the `i`-th spatial dimension of the original image.
"""
function Images.component_centroids(c::CC, coords)
    n = length(c)
    D = length(coords)
    s = ntuple(_ -> 0.0, D)
    for p in c
        s = s .+ ntuple(i -> coords[i][p[i]], D)
    end
    return s ./ n
end