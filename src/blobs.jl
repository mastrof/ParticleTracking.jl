export AbstractBlob, Blob, DummyBlob

"""
    AbstractBlob{T,S,N}
Abstract interface for blob types
"""
abstract type AbstractBlob{T,S,N} end
function Base.show(io::IO, blob::AbstractBlob)
    l = location(blob)
    σ = scale(blob)
    a = amplitude(blob)
    m0 = zeroth_moment(blob)
    m2 = second_moment(blob)
    print(io, 
        "Blob(location=$l, σ=$σ, amplitude=$a, m₀=$m0, m₂=$m2)"
    )
end


"""
    struct BlobRaw{T,S,N} <: AbstractBlob{T,S,N}
        location::CartesianIndex{N}
        σ::S
        amplitude::T
        m0::T
        m2::T
        intensity_map
    end
"""
struct BlobRaw{T,S,N} <: AbstractBlob{T,S,N}
    location::CartesianIndex{N}
    σ::S
    amplitude::T
    m0::T
    m2::T
    intensity_map
end

function BlobRaw(
    location::CartesianIndex{N}, σ::S, amplitude::T, img::AbstractArray{T,N}
) where {T,S,N}
    img_pad, m0, m2 = image_moments(location, σ, img)
    BlobRaw(location, σ, amplitude, m0, m2, img_pad)
end

"""
    struct Blob{T,S,N} <: AbstractBlob{T,S,N}
        location::NTuple{N,Float64}
        location_raw::CartesianIndex{N}
        σ::S
        amplitude::T
        m0::T
        m2::T
        intensity_map
    end
"""
struct Blob{T,S,N} <: AbstractBlob{T,S,N}
    location::NTuple{N,Float64}
    location_raw::CartesianIndex{N}
    σ::S
    amplitude::T
    m0::T
    m2::T
    intensity_map
end

function Blob(blob::AbstractBlob{T,S,N}, offsets::NTuple{N,Float64}) where {T,S,N}
    raw_location = location(blob)
    Blob(Tuple(raw_location) .+ offsets, raw_location, scale(blob),
        amplitude(blob), zeroth_moment(blob), second_moment(blob), intensity_map(blob)
    )
end


function DummyBlob(BlobType::Type{Blob{T,S,N}}) where {T,S,N}
    BlobType(
        ntuple(_ -> T(NaN), N),
        CartesianIndex(ntuple(_ -> 0, N)),
        S(0 for _ in 1:N),
        zero(T),
        zero(T),
        zero(T),
        []
    )
end
