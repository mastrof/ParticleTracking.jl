module ParticleTracking

using JLD2
using LinearAlgebra, StatsBase, Statistics
using OffsetArrays
using Hungarian
using Images

# shorthand notations
CI = CartesianIndex
CC{N} = Vector{CI{N}} # connected component
Stack{T,N} = Vector{Array{T,N}}
Stack2{T} = Stack{T,2}
Stack3{T} = Stack{T,3}

end
