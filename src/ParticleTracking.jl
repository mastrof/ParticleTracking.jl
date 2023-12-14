module ParticleTracking

using JLD2
using LinearAlgebra, StatsBase, Statistics
using OffsetArrays
using Hungarian
using Images
using BlobDetection

export Trajectory

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


include("linking.jl")
include("trajectory_api.jl")

end
