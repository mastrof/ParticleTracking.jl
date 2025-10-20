export track_blobs

"""
    track_blobs(blobs)

Track `blobs` over time to produce a vector of trajectories by solving the
associated minimum-cost flow problem.

**Keyword arguments**
- `minlife::Integer`: trajectories shorter than this value are discarded; defaults to 1
- `memory::Integer`: number of frames over which a blob is allowed to disappear; defaults to 0
- `cost::Function`: cost function used to link blobs; must accept two
  `AbstractBlob`s as arguments (along with any desired kwarg) and return a non-negative
  real number.
- `maxdist::Real`: maximum distance that a blob is allowed to travel in a single frame; defaults to Inf
- `gapscale::Function`: scaling function to reweight allowed distance over detection gaps, i.e.,
  the maxdist over n frames is maxdist * gapscale(n); defaults to sqrt.
- `k::Integer`: lookback window for velocity prediction; defaults to 0 (no prediction)
- `w::Real`: weight (must be within [0,1]) assigned to velocity prediction wrt spatial linking;
  defaults to 0.

Other arbitrary kwargs can be provided, they will be passed directly to the `cost` function.
"""
function track_blobs(
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    minlife::Integer = 1,
    memory::Integer = 1,
    cost::Function = QuadraticCost,
    maxdist::Real = Inf,
    gapscale::Function = sqrt,
    k::Integer = 0,
    w::Real = 0.0,
    kwargs...
)
    @assert minlife >= 1
    @assert memory >= 0
    @assert maxdist > 0
    @assert k >= 0
    @assert 0 <= w <= 1
    links = tracking(blobs, memory, cost, maxdist, gapscale, k, w; kwargs...)
    trajectories = build_trajectories(links, blobs; minlife)
end

function tracking(
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob{T,N}}},
    memory::Integer,
    cost::Function,
    maxdist::Real,
    gapscale::Function,
    k::Integer,
    w::Real;
    kwargs...
) where {T,N}
    num_frames = length(blobs)
    linked_from = [falses(length(b)) for b in blobs]
    linked_to = [falses(length(b)) for b in blobs]
    all_links = DataFrame(from_t=Int[], from_idx=Int[], to_t=Int[], to_idx=Int[])
    # keep track of parent-child relationships for velocity computations
    ancestor_map = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    # precompute max allowed costs
    BT = eltype(eltype(blobs))
    maxcost = [
        isinf(maxdist) ? maxdist : evaluate_maxcost(BT, maxdist, gapscale(r), cost; kwargs...)
        for r in 1:(1+memory)
    ]
    # using outer r loop and inner t loop ensures that
    # links with smaller frame gap are given priority
    for r in 1:(1+memory)
        for t in 1:num_frames-r
            #== velocity calculation ==#
            velocities = Dict{Tuple{Int,Int}, NTuple{N,Float64}}()
            estimate_velocities!(velocities, ancestor_map, blobs, t, k)
            predictive_cost = PredictiveCost(cost, velocities, float(w))
            #== linking ==#
            # find link candidates between frames t and t+r
            from_indices = findall(.!linked_from[t])
            to_indices = findall(.!linked_to[t+r])
            if isempty(from_indices) || isempty(to_indices)
                continue
            end
            from_blobs = blobs[t][from_indices]
            to_blobs = blobs[t+r][to_indices]
            # setup linking problem
            Nt = length(from_blobs)
            Nr = length(to_blobs)
            linking_matrix = OffsetMatrix(zeros(Int, Nt+1, Nr+1), 0:Nt, 0:Nr)
            cost_matrix = OffsetMatrix(zeros(Float64, Nt+1, Nr+1), 0:Nt, 0:Nr)
            # populate cost matrix
            evaluate_costs!(
                cost_matrix, from_blobs, to_blobs, r,
                predictive_cost, maxcost[r], t, from_indices;
                kwargs...
            )
            # linearity of cost functional ensures existence of unique optimum
            # but good initialization of links can improve performance
            initialize_links!(linking_matrix, cost_matrix, from_blobs, to_blobs, t, r)
            # rearrange links until optimal configuration is found
            optimize_links!(linking_matrix, cost_matrix)
            # establish connections from linking matrix
            new_links = makelinks(linking_matrix, from_indices, to_indices, t, r)
            for link in new_links
                ancestor_map[(link.to_t, link.to_idx)] = (link.from_t, link.from_idx)
                linked_from[t][link.from_idx] = true
                linked_to[t+r][link.to_idx] = true
                push!(all_links, link)
            end
        end
    end
    all_links
end

# construct trajectory objects from dataframe of links
function build_trajectories(
    links_df::DataFrame,
    blobs::AbstractVector{<:AbstractVector{B}};
    minlife::Integer=1
) where {B<:AbstractBlob}
    # lookup map for traversal
    links_map = Dict{Tuple{Int,Int}, Tuple{Int, Int}}()
    to_nodes = Set{Tuple{Int,Int}}()
    for row in eachrow(links_df)
        from_node = (row.from_t, row.from_idx)
        to_node = (row.to_t, row.to_idx)
        links_map[from_node] = to_node
        push!(to_nodes, to_node)
    end
    # find starting points
    from_nodes = Set(keys(links_map))
    start_nodes = setdiff(from_nodes, to_nodes)
    trajectories = Trajectory{Int, B}[]
    visited_nodes = Set{Tuple{Int, Int}}()
    # build trajectory
    for start_node in start_nodes
        (start_node in visited_nodes) && continue
        current_times = Int[]
        current_blobs = B[]
        current_node = start_node
        # walk the path until end
        while true
            t, idx = current_node
            push!(current_times, t)
            push!(current_blobs, blobs[t][idx])
            push!(visited_nodes, current_node)
            # find next node
            next_node = get(links_map, current_node, nothing)
            isnothing(next_node) && break
            current_node = next_node
        end
        trajectory = Trajectory(current_times, current_blobs)
        if length(trajectory) >= minlife
            push!(trajectories, trajectory)
        end
    end
    trajectories
end
