module ParticleTracking

using OffsetArrays
using DataFrames
using ImageFiltering
using PaddedViews
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
using LinearAlgebra: norm, dot
include("linking_cost.jl")
include("prediction.jl")
include("linking.jl")
include("tracking.jl")

#== trackpy ==#
include("trackpy.jl")

#== Visualization ==#
export explore_blobs, explore_tracking

"""
    explore_blobs(video; kwargs...)
Open a window with an interactive preview of the blob detection on `video`.
The UI allows to move through all the frames of the `video` while tuning
the blob scale and amplitude threshold for detection.
**Only available upon loading GLMakie**.
"""
function explore_blobs end

"""
    explore_tracking(video, blobs; kwargs...)
Open a window with an interactive preview of the tracking on `video` based
on the detected `blobs`.
"""
function explore_tracking end

end
