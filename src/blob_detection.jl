export detect_blobs

"""
    detect_blobs(img, σscales; edges::Tuple=(true, false, ...), σshape::Tuple=(1, ...), rthresh=0.001) -> Vector{Blob}
    detect_blobs(stack, σscales; ...)

Find "blobs" in a N-dimensional image using the negative Laplacian of Gaussian (LoG) filter
for the specified set of scales (`σscales`).
If a stack (i.e. a vector of N-dimensional images) is passed instead,
blobs are identified in each frame individually.
The algorithm searches for places where the filtered image (for a particular σ) is at a peak
compared to all spatially and σ-adjacent voxels, where σ is `σscales[i] * σshape` for some `i`.
By default, `σshape` is an ntuple of 1s.

The optional `edges` argument controls whether peaks on the edges are included.
`edges` can be `true` or `false`, or a N+1-tuple in which the first entry controls whether
edge-σ values are eligible to serve as peaks, and the remaining N control each of the N
dimensions of `img`.

`rthresh` controls the minimum amplitude of peaks in the -LoG-filtered image, as a fraction
of `maximum(abs, img)`.

This method is based on a minor modification of the `blob_LoG` algorithm implemented in
the ImageFiltering package.
"""
function detect_blobs end
function detect_blobs(stack::AbstractVector{<:AbstractArray{T,N}}, σscales;
    kwargs...
) where {T<:Real,N}
    map(img -> detect_blobs(img, σscales; kwargs...), stack)
end
function detect_blobs(img::AbstractArray{T,N}, σscales;
    edges::Union{Bool,Tuple{Bool,Vararg{Bool,N}}} = (true, ntuple(d -> false, Val(N))...),
    σshape::NTuple{N,Real} = ntuple(d -> 1, Val(N)),
    rthresh::Real = 1 // 1000
) where {T<:Real,N}
    if edges isa Bool
        edges = (edges, ntuple(d -> edges, Val(N))...)
    end
    sigmas = sort(σscales)
    # apply -LoG filter for each scale
    # and cap all negative values to 0
    img_LoG = max.(multiLoG(img, sigmas, σshape), zero(T))
    colons = ntuple(d -> Colon(), Val(N))
    maxima = findlocalmaxima(img_LoG; edges)
    if !iszero(rthresh)
        imgmax = maximum(abs, img)
        [refine(BlobRaw(
            CartesianIndex(Base.tail(x.I)),
            sigmas[x[1]].*σshape,
            img_LoG[x],
            img_LoG[x[1], colons...]
            )) for x in maxima if img_LoG[x] > rthresh*imgmax
        ]
    else
        [refine(BlobRaw(
            CartesianIndex(Base.tail(x.I)),
            sigmas[x[1]].*σshape,
            img_LoG[x],
            img_LoG[x[1], colons...]
            )) for x in maxima
        ]
    end
end

function multiLoG(img::AbstractArray{T,N}, sigmas, σshape) where {T,N}
    issorted(sigmas) || error("sigmas must be sorted")
    img_LoG = similar(img, float(eltype(T)), (Base.OneTo(length(sigmas)), axes(img)...))
    colons = ntuple(d -> Colon(), Val(N))
    @inbounds for (isigma, σ) in enumerate(sigmas)
        LoG_slice = @view img_LoG[isigma, colons...]
        imfilter!(LoG_slice, img, Kernel.LoG(ntuple(i -> σ*σshape[i], Val(N))), "reflect")
        LoG_slice .*= -σ^2
    end
    return img_LoG
end
