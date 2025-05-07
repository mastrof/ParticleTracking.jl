export Trajectory

struct Trajectory{T<:Integer,B<:AbstractBlob}
    times::Vector{T}
    blobs::Vector{B}
end
function Base.show(io::IO, trajectory::Trajectory)
    l = length(trajectory)
    x0 = startpoint(trajectory)
    x1 = endpoint(trajectory)
    t0, t1 = extrema(timestamps(trajectory))
    print(io,
        "Trajectory($l points from [x=$x0, t=$t0] to [x=$x1, t=$t1]"
    )
end

#== build trajectory from blobs and links ==#
function Trajectory(
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}},
    connections::AbstractVector{<:Integer}
)
    times = findall(!iszero, connections)
    trajectory_blobs = [blobs[t][connections[t]] for t in times]
    Trajectory(times, trajectory_blobs)
end

#== Iteration and indexing ==#
Base.iterate(X::Trajectory, state=1) = iterate(X.blobs, state)
function Base.getindex(X::Trajectory, i::Integer)
    X.blobs[i]
end
function Base.getindex(X::Trajectory, inds)
    [X[i] for i in inds]
end
# call trajectory as a function
# works like indexing if a detection exists at that timepoint
# otherwise returns a dummy blob with NaN location
function (traj::Trajectory{T,B})(i::Integer) where {T,S,N,B<:AbstractBlob{S,N}}
    ts = timestamps(traj)
    t1 = ts[1]
    t2 = ts[end]
    t = findfirst(ts .== i)
    if isnothing(t)
        return B(; location=ntuple(_ -> S(NaN), N))
    else
        return traj.blobs[t]
    end
end
