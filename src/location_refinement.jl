export BlobRefined, refine, location_refined

struct BlobRefined{T,S,N} <: AbstractBlob{T,S,N}
    location::NTuple{N,T}
    location_raw::CartesianIndex{N}
    σ::S
    amplitude::T
    m0::T
    m2::T
    intensity_map
end
location_raw(blob::BlobRefined) = blob.location_raw

function BlobRefined(blob::AbstractBlob{T,S,N}, offsets::NTuple{N,T}) where {T,S,N}
    raw_location = location(blob)
    BlobRefined(Tuple(raw_location) .+ offsets, raw_location, scale(blob),
        amplitude(blob), zeroth_moment(blob), second_moment(blob), intensity_map(blob)
    )
end

"""
    refine(blob::AbstractBlob)
"""
function refine(blob::AbstractBlob)
    ε = offsets(blob)
    return BlobRefined(blob, ε)
end
refine(blob::BlobRefined) = blob # do nothing if position is already refined

function offsets(blob::AbstractBlob{T,S,2}) where {T,S}
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
