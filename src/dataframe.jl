export to_df

"""
    to_df(tracks::AbstractVector{<:Trajectory})
Convert a set of trajectories into a dataframe.
"""
function to_df(tracks::AbstractVector{<:Trajectory})
    vcat(map(to_df, tracks, eachindex(tracks))...)
end
function to_df(track::Trajectory{S,<:AbstractBlob{T,2}}, id::Integer) where {S,T}
    frame = timestamps(track)
    particle = fill(id, length(frame))
    pos = location(track)
    x = first.(pos)
    y = last.(pos)
    m0 = zeroth_moment(track)
    m2 = second_moment(track)
    s = first.(scale(track))
    DataFrame(; frame, particle, x, y, mass=m0, m2, size=s)
end
function to_df(track::Trajectory{S,<:AbstractBlob{T,3}}, id::Integer) where {S,T}
    frame = timestamps(track)
    particle = fill(id, length(frame))
    pos = location(track)
    x = first.(pos)
    y = getindex.(pos, 2)
    z = last.(pos)
    m0 = zeroth_moment(track)
    m2 = second_moment(track)
    s = first.(scale(track))
    DataFrame(; frame, particle, x, y, z, mass=m0, m2, size=s)
end

"""
    to_df(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}})
Convert a vector of blob detections into a dataframe.
"""
function to_df(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}})
    vcat(map(to_df, blobs, eachindex(blobs))...)
end
function to_df(blobs::AbstractVector{<:AbstractBlob{T,2}}, t::Integer) where T
    frame = fill(t, length(blobs))
    pos = location.(blobs)
    x = first.(pos)
    y = last.(pos)
    m0 = zeroth_moment.(blobs)
    m2 = second_moment.(blobs)
    s = first.(scale.(blobs))
    DataFrame(; frame, x, y, mass=m0, m2, size=s)
end
function to_df(blobs::AbstractVector{<:AbstractBlob{T,3}}, t::Integer) where T
    frame = fill(t, length(blobs))
    pos = location.(blobs)
    x = first.(pos)
    y = getindex.(pos, 2)
    z = last.(pos)
    m0 = zeroth_moment.(blobs)
    m2 = second_moment.(blobs)
    s = first.(scale.(blobs))
    DataFrame(; frame, x, y, z, mass=m0, m2, size=s)
end
