module ParticleTracking

using OffsetArrays
using ImageFiltering
using PaddedViews
using MakieCore
using GeometryBasics


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
using Requires
function __init__()
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("gui_blobs.jl")
end

end
