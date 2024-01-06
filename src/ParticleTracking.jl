module ParticleTracking

using OffsetArrays
using ImageFiltering
using PaddedViews
using MakieCore, GeometryBasics


#== Type interface ==#
include("blobs.jl")
include("blob_api.jl")
include("trajectories.jl")
include("trajectory_api.jl")

#== Blob detection ==#
include("image_moments.jl")
include("blob_detection.jl")
include("location_refinement.jl")

#== Tracking ==#
#include("linking.jl")
using LinearAlgebra: norm, dot
include("linking_cost.jl")
include("flow_linking.jl")

#== Visualization ==#
include("makie_recipes.jl")

end
