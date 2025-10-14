# draw blobs as scatter points
function Makie.convert_arguments(
    P::Type{<:Scatter},
    x::AbstractVector{<:AbstractBlob}
)
    (Point2f.(location.(x)),)
end

# draw blobs as circles
function Makie.convert_arguments(
    P::Type{<:Poly},
    x::AbstractVector{<:AbstractBlob},
    r::Union{<:Real,<:AbstractVector{<:AbstractBlob}}=1
)
    (Circle.(Point2f.(location.(x)), r),)
end

# draw a sequence of blobs as a trajectory
function Makie.convert_arguments(
    P::Type{<:Lines},
    x::AbstractVector{<:AbstractBlob}
)
    (Point2f.(location.(x)),)
end

# draw a trajectory up to time t with a tail of length t0
function Makie.convert_arguments(
    P::Type{<:Lines},
    x::Trajectory,
    t::Integer,
    t0::Integer
)
    (Point2f.(location.(x.(max(1,t-t0):t))),)
end
# draw entire trajectory
Makie.convert_arguments(P::Type{<:Lines}, x::Trajectory) =
    (Point2f.(location(x)),)
