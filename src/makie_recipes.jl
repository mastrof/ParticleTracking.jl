# draw a trajectory up to time t with a tail of length t0
MakieCore.convert_arguments(P::Type{<:Lines}, x::Trajectory, t::Integer, t0::Integer) =
    (Point2f.(location.(x[max(1,t-t0):t])),)
# draw entire trajectory
MakieCore.convert_arguments(P::Type{<:Lines}, x::Trajectory) =
    (Point2f.(location(x)),)

# draw a sequence of blobs as a trajectory
MakieCore.convert_arguments(P::Type{<:Lines}, x::AbstractVector{<:AbstractBlob}) =
    (Point2f.(location.(x)),)

# draw blobs as scatter points
MakieCore.convert_arguments(P::Type{<:Scatter}, x::AbstractVector{<:AbstractBlob}) =
    (Point2f.(location.(x)),)
