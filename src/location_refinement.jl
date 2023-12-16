export refine, offsets

"""
    refine(blob::AbstractBlob)
"""
function refine(blob::AbstractBlob)
    ε = offsets(blob)
    return BlobRefined(blob, ε)
end
refine(blob::BlobRefined) = blob # do nothing if position is already refined

function offsets(blob::Blob{T,S,2}) where {T,S}
    x, y = location(blob).I
    σx, σy = scale(blob)
    m0 = zeroth_moment(blob)
    I = intensity_map(blob)
    εx = zero(float(T))
    εy = zero(float(T))
    for j in axes(I,2), i in axes(I,1)
        dy = y-j
        dx = x-i
        εx += dx * I[x+dx, y+dy]
        εy += dy * I[x+dx, y+dy]
    end
    εx /= m0
    εy /= m0
    if abs(εx) > 0.5 || abs(εy) > 0.5
        new_location = CartesianIndex(x + 1*sign(εx), y + 1*sign(εy))
        new_blob = Blob(new_location, scale(blob), amplitude(blob),
            zeroth_moment(blob), second_moment(blob), intensity_map(blob)
        )
        offsets(new_blob)
    end
    return (εx, εy)
end
offsets(blob::BlobRefined{T,S,N}) where {T,S,N} = ntuple(_ -> 0.0, N)
