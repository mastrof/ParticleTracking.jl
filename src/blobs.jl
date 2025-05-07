export AbstractBlob, Blob, DummyBlob

"""
    AbstractBlob{T,S,N}
Abstract interface for blob types
"""
abstract type AbstractBlob{T,N} end
function Base.show(io::IO, blob::AbstractBlob)
    l = location(blob)
    print(io, "Blob(location=$l)")
end


"""
    struct BlobRaw{T,N} <: AbstractBlob{T,N}
        location::CartesianIndex{N}
        σ::NTuple{N,<:Number}
        amplitude::T
        m0::T
        m2::T
        intensity_map
    end
"""
@kwdef struct BlobRaw{T,N} <: AbstractBlob{T,N}
    location::CartesianIndex{N}
    σ::NTuple{N,<:Number} = ntuple(zero, N)
    amplitude::T = zero(T)
    m0::T = zero(T)
    m2::T = zero(T)
    intensity_map = Array{T,N}(undef, ntuple(zero, N)...)
end

function BlobRaw(
    location::CartesianIndex{N},
    σ::NTuple{N,<:Number},
    amplitude::T,
    img::AbstractArray{T,N}
) where {T,N}
    img_pad, m0, m2 = image_moments(location, σ, img)
    BlobRaw(location, σ, amplitude, m0, m2, img_pad)
end

"""
    struct Blob{T,N} <: AbstractBlob{T,N}
        location::NTuple{N,Float64}
        location_raw::CartesianIndex{N}
        σ::NTuple{N,<:Number}
        amplitude::T
        m0::T
        m2::T
        intensity_map
    end
"""
@kwdef struct Blob{T,N} <: AbstractBlob{T,N}
    location::NTuple{N,Float64}
    location_raw::CartesianIndex{N} = CartesianIndex(
        any(isnan, location) ? ntuple(zero, N) : round.(Int, location)
    )
    σ::NTuple{N,<:Number} = ntuple(zero, N)
    amplitude::T = zero(T)
    m0::T = zero(T)
    m2::T = zero(T)
    intensity_map = Array{T,N}(undef, ntuple(zero, N)...)
end

function Blob(
    blob::AbstractBlob{T,N},
    offsets::NTuple{N,<:Number}
) where {T,N}
    raw_location = location(blob)
    Blob(Tuple(raw_location) .+ offsets, raw_location, scale(blob),
        amplitude(blob), zeroth_moment(blob), second_moment(blob), intensity_map(blob)
    )
end
