export refine, offsets

"""
    refine(blob::AbstractBlob)
Estimate the real location of `blob` through subpixel localization.
"""
function refine(blob::AbstractBlob)
    ε = offsets(blob)
    return Blob(blob, ε)
end
refine(blob::Blob) = blob # do nothing if position is already refined

function offsets(blob::BlobRaw{T,S,N}) where {T,S,N}
    m0 = zeroth_moment(blob)
    I = intensity_map(blob)
    x = location(blob).I
    ε = zeros(N)
    for c in CartesianIndices(I)
        delta = x .- c.I
        ε .+= delta .* I[CartesianIndex(x .+ delta)]
    end
    ε ./= m0
    # BUG: was bugged with previous 2D implementation; should work now but needs check.
    # When new_location is set it does not coincide with the center of intensity map;
    # on the following iteration indexing is messed up and goes out of bounds
    #==
    if any(abs(d) > 0.5 for d in ε)
        new_location = CartesianIndex(@. x + Int(sign(ε)))
        new_blob = BlobRaw(new_location, scale(blob), amplitude(blob),
            zeroth_moment(blob), second_moment(blob), intensity_map(blob)
        )
        offsets(new_blob)
    end
    ==#
    return Tuple(ε)
end
offsets(blob::Blob{T,S,N}) where {T,S,N} = ntuple(_ -> 0.0, N)
