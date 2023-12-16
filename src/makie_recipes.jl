function MakieCore.convert_arguments(P::Type{<:Lines}, x::Trajectory, t::Integer, t0::Integer)
    (Point2f.(location.(x[max(1,t-t0):t])),)
end
MakieCore.convert_arguments(P::Type{<:Lines}, x::AbstractVector{<:AbstractBlob}) =
    (Point2f.(location.(x)),)

MakieCore.convert_arguments(P::Type{<:Scatter}, x::AbstractVector{<:AbstractBlob}) =
    (Point2f.(location.(x)),)
