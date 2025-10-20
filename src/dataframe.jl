export to_df

"""
    to_df(tracks::AbstractVector{<:Trajectory})
Convert a set of trajectories into a dataframe.
"""
function to_df(tracks::AbstractVector{<:Trajectory})
    df = DataFrame(;
        frame=Int[],
        particle=Int[],
        x=Float64[],
        y=Float64[],
        m0=Float64[],
        m2=Float64[],
        size=Float64[],
    )
    foreach((track, id) -> append!(df, to_df(track, id)), tracks, 1:length(tracks))
    df
end
"""
    to_df(track, id)
Convert a single track to a dataframe and assign it the given id.
"""
function to_df(track::Trajectory, id::Integer)
    frame = timestamps(track)
    particle = fill(id, length(frame))
    pos = location(track)
    x = first.(pos)
    y = last.(pos)
    m0 = zeroth_moment(track)
    m2 = second_moment(track)
    s = first.(scale(track))
    DataFrame(; frame, particle, x, y, m0, m2, size=s)
end
