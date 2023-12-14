module ParticleTracking

using JLD2
using LinearAlgebra, StatsBase, Statistics
using OffsetArrays
using Hungarian
using Images
using ImageFiltering
using PaddedViews
using ColorTypes

export AbstractBlob, Trajectory

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

struct Trajectory{T<:Integer,B<:AbstractBlob}
    times::Vector{T}
    blobs::Vector{B}
end

function Trajectory(
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}},
    connections::AbstractVector{<:Integer}
)
    times = findall(!iszero, connections)
    trajectory_blobs = [blobs[t][connections[t]] for t in times]
    Trajectory(times, trajectory_blobs)
end


#== Blob detection ==#
include("image_moments.jl")
include("blob_api.jl")
include("blobs.jl")
include("location_refinement.jl")

#== Tracking ==#
include("linking.jl")
include("trajectory_api.jl")

end
