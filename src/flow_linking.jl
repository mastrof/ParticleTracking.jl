export blobtracking

"""
    blobtracking(blobs; minlife=1, memory=1, cost=QuadraticCost, maxdist=Inf, kwargs...)

Track `blobs` over time to produce a vector of trajectories by solving the
associated minimum-cost flow problem.

**Keyword arguments**
- `minlife::Integer`: trajectories shorter than this value are discarded
- `memory::Integer`: number of frames over which a blob is allowed to disappear
- `cost::Function`: cost function used to link blobs; must accept two
  `AbstractBlob`s as arguments (along with any desired kwarg) and return a non-negative
  real number.
- `maxdist::Real`: maximum distance that a blob is allowed to travel in a single frame

Other arbitrary kwargs can be provided, they will be passed directly to the `cost` function.
"""
function blobtracking(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    minlife::Integer = 1,
    memory::Integer = 1,
    cost::Function = QuadraticCost,
    maxdist::Real = Inf,
    kwargs...
)
    @assert minlife >= 1
    @assert memory >= 1
    @assert maxdist > 0
    links = mincostflow(blobs; memory, cost, maxdist, kwargs...)
    trajectories = buildtraj(blobs, links; minlife)
end

function mincostflow(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    memory::Integer, cost::Function, maxdist::Real, kwargs...
)
    times = eachindex(blobs)
    T = eltype(eltype(blobs))
    # precompute max allowed costs for linking
    maxcost = [
        isinf(maxdist) ? maxdist : evaluate_maxcost(T, maxdist, r, cost; kwargs...)
        for r in 1:memory
    ]
    [mincostflow(blobs, t, memory; cost, maxcost, kwargs...) for t in times]
end

function mincostflow(blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}},
    t::Integer, memory::Integer;
    cost::Function, maxcost::AbstractVector{<:Real}, kwargs...
)
    tmax = length(blobs)
    from = blobs[t]
    Nt = length(from)
    to = blobs[t+1:min(t+memory,tmax)]
    Nr = length.(to)
    # row and column 0 are introduced for blob creation/annihilation
    linking_matrix = [OffsetMatrix(zeros(Int, Nt+1, n+1), 0:Nt, 0:n) for n in Nr]
    cost_matrix = [OffsetMatrix(zeros(Nt+1, n+1), 0:Nt, 0:n) for n in Nr]
    # linearity of cost functional ensures existence of unique optimum
    # but a good initialization of the links can improve performance
    initialize_links!(linking_matrix, cost_matrix, from, to, t; cost, maxcost, kwargs...)
    # rearrange links until optimal configuration is found
    optimize_links!(linking_matrix, cost_matrix, from, to)
    # transform linking matrix into a sequence of connections
    links = makelinks(linking_matrix, from, to, t)
end

function initialize_links!(G::AbstractVector{TG}, C::AbstractVector{TC},
    from::AbstractVector{<:AbstractBlob}, to::AbstractVector{VB}, t::Integer;
    cost::Function, maxcost::AbstractVector{<:Real}, kwargs...
) where {TG<:OffsetMatrix,TC<:OffsetMatrix,VB<:AbstractVector{<:AbstractBlob}}
    R = eachindex(to)
    for r in R
        evaluate_costs!(C[r], from, to[r], r, cost, maxcost[r]; kwargs...)
        initialize_links!(G[r], C[r], from, to[r], t, r)
    end
end

function initialize_links!(G::OffsetMatrix, C::OffsetMatrix,
    from::AbstractVector{<:AbstractBlob}, to::AbstractVector{<:AbstractBlob},
    t::Integer, r::Integer
)
    js_occupied = zero(axes(G, 2))
    # initial guesses for min cost links
    for i in axes(G,1)[1:end]
        # TODO: add check to avoid throwing errors when empty
        c, j = findmin(C[i,j] for j in axes(G,2)[1:end])
        if js_occupied[j] == 0 && !isinf(C[i,j])
            G[i,j] = 1
            if j != 0
                js_occupied[j] = 1
            end
        else
            G[i,0] = 1
        end
    end
    # link all unassigned blobs to dummies
    for j in axes(G,2)[1:end]
        if sum(G[:,j]) == 0
            G[0,j] = 1
        end
    end
end

function optimize_links!(
    G::AbstractVector{<:OffsetMatrix}, C::AbstractVector{<:OffsetMatrix},
    from::AbstractVector{<:AbstractBlob},
    to::AbstractVector{<:AbstractVector{<:AbstractBlob}},
)
    for r in eachindex(to)
        Z = zero(C[r]) .- eltype(C[r])(Inf) # preallocate reduced cost matrix
        while any(Z .< 0)
            optimize_links!(G[r], Z, C[r], from, to[r])
        end
    end
end
function optimize_links!(
    G::OffsetMatrix, Z::OffsetMatrix, C::OffsetMatrix,
    from::AbstractVector{<:AbstractBlob}, to::AbstractVector{<:AbstractBlob},
)
    Z .= zero(eltype(Z)) # reset all reduced costs
    # look for non-linked pairs with finite cost
    indices = findall(idx -> G[idx] == 0 && !isinf(C[idx]), CartesianIndices(G))
    minZ = 0.0
    minIJKL = [0, 0, 0, 0]
    for (I,J) in Tuple.(indices)
        if I > 0 && J > 0 # swap
            L = findfirst(G[I,j] == 1 for j in axes(G,2))
            K = findfirst(G[i,J] == 1 for i in axes(G,1))
            Z[I,J] = C[I,J] - C[I,L] - C[K,J] + C[K,L]
        elseif I == 0 && J > 0 # creation
            L = 0
            K = findfirst(G[i,J] == 1 && i > 0 for i in axes(G,1))
            Z[0,J] = C[0,J] - C[K,J] + C[K,0]
        elseif I > 0 && J == 0 # annihilation
            L = findfirst(G[I,j] == 1 && j > 0 for j in axes(G,2))
            K = 0
            Z[I,0] = C[I,0] - C[I,L] + C[0,L]
        end
        if Z[I,J] < minZ
            minZ = Z[I,J]
            minIJKL .= (I,J,K,L)
        end
    end
    # update links with minimal reduced cost if any
    if minZ < 0
        I,J,K,L = minIJKL
        G[I,J] = 1
        if I != 0
            G[I,L] = 0
        end
        if J != 0
            G[K,J] = 0
        end
        G[K,L] = 1
    end
    nothing
end

function makelinks(
    G::AbstractVector{<:OffsetMatrix},
    from::AbstractVector{<:AbstractBlob},
    to::AbstractVector{<:AbstractVector{<:AbstractBlob}},
    t::Integer
)
    links = Dict{Int, Tuple{Int, Int}}() # i => (j, t+r)
    is_disconnected = [true for _ in from] # no links connected yet
    # check connections in order through successive frames
    for r in eachindex(to)
        connections = Tuple.(findall(G[r] .== 1))
        for (i,j) in connections
            if i>0 && is_disconnected[i] && j!=0
                links[i] = (j, t+r)
                is_disconnected[i] = false
            end
        end
        !any(is_disconnected) && break
    end
    return links
end

function buildtraj(
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}},
    links::AbstractVector{<:AbstractDict};
    minlife=1
)
    connections = joinlinks(links; minlife)
    trajectories = [buildtraj(blobs, c) for c in connections]
end

function buildtraj(blobs, connections)
    # each element of connections is (i,t)
    traj_blobs = [blobs[t][i] for (i,t) in connections]
    traj_times = last.(connections)
    Trajectory(traj_times, traj_blobs)
end

function joinlinks(linkmap::AbstractVector{<:AbstractDict}; minlife=1)
    linkmap_copy = deepcopy(linkmap)
    tracks = Vector{Tuple{Int,Int}}[]
    for t in eachindex(linkmap)
        for i in keys(linkmap_copy[t])
            traj = joinlinks(linkmap_copy, t, i)
            push!(tracks, traj)
            for (k,s) in traj
                s <= length(linkmap) && delete!(linkmap_copy[s], k)
            end
        end
    end
    filter(s -> length(s) >= minlife, tracks)
end

function joinlinks(linkmap, t, i)
    tmax = length(linkmap) - t
    j = i
    traj = Tuple{Int,Int}[(i,t)]
    sizehint!(traj, tmax)
    s = t
    while s < tmax
        k, s = traj[end]
        if s <= tmax && haskey(linkmap[s], k)
            nextid, nextframe = linkmap[s][k]
            push!(traj, (nextid, nextframe))
        else
            break
        end
    end
    traj
end
