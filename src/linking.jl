export track, link

function track(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    cost::Function = Quadratic,
    minlife::Int = 1,
    kwargs...
)
    links = link(blobs; cost, kwargs...)
    connections = connect(links; minlife)
    trajectories = [Trajectory(blobs, c) for c in connections]
end

function link(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    cost::Function = Quadratic,
    kwargs...
)
    indices = eachindex(blobs)[2:end]
    [link(blobs[t-1], blobs[t]; cost, kwargs...) for t in indices]
end

function link(
    blobs_past::AbstractVector{<:AbstractBlob},
    blobs_next::AbstractVector{<:AbstractBlob};
    cost::Function = Quadratic,
    kwargs...
)
    weights = [cost(A, B; kwargs...) for A in blobs_past, B in blobs_next]
    assignment, _ = hungarian(weights)
    Dict(i => j for (i,j) in enumerate(assignment))
end

function connect(links; minlife=1)
    links_copy = deepcopy(links)
    tracks = Vector{Int}[]
    for t in eachindex(links)
        for i in keys(links_copy[t])
            traj = connect(links_copy, t, i)
            push!(tracks, traj)
            for (s,k) in enumerate(traj[1:end-1])
                delete!(links_copy[s], k)
            end
        end
    end
    filter(s -> count(s .> 0) >= minlife, tracks)
end

function connect(links, t, i)
    tmax = length(links) - t
    j = i
    traj = [repeat([0], t-1); Int[i]]
    s = -1
    while s < tmax
        k = traj[end]
        s += 1
        j = links[t+s][k]
        j == 0 && break
        push!(traj, j)
    end
    vcat(traj, repeat([0], length(links)-length(traj)+1))
end

function Quadratic(A::AbstractBlob, B::AbstractBlob; 
    maxdist::Real = Inf,
    g0::Real = 1.0,
    g2::Real = 1.0
)
    lA = Tuple(location(A))
    lB = Tuple(location(B))
    m0A = zeroth_moment(A)
    m0B = zeroth_moment(B)
    m2A = second_moment(A)
    m2B = second_moment(B)
    s1 = sum(abs2.(lA .- lB))
    if s1 > maxdist^2
        return missing
    end
    s2 = g0*(m0A - m0B)^2
    s3 = g2*(m2A - m2B)^2
    return s1 + s2 + s3
end
