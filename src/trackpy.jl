export trackpydf_to_blobs

"""
    trackpydf_to_blobs(df::AbstractDataFrame)
Convert a dataframe containing particle detections from trackpy
to a collection of blobs.
Returns a vector of vector of blobs
(a vector of blobs for each frame).

Currently only works for 2-dimensional systems.
"""
function trackpydf_to_blobs(df::AbstractDataFrame)
    N = 2 # TODO: allow 3d detections
    T = Float64
    gdf = groupby(sort(df, :frame), :frame)
    [
        map(
            row -> Blob{T,N}(;
                location=(row.x, row.y),
                σ=(row.size, row.size),
                amplitude=row.signal,
                m0=row.mass
            ),
            eachrow(g)
        )
        for g in gdf
    ]
end
