isnantuple(t::Tuple) = any(isnan, t)

function ParticleTracking.explore_tracking(
    video::AbstractArray{T,3},
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    kwargs...
) where T
    explore_tracking(eachslice(video; dims=3), blobs; kwargs...)
end
function ParticleTracking.explore_tracking(
    video::AbstractVector{<:AbstractMatrix{T}},
    blobs::AbstractVector{<:AbstractVector{<:AbstractBlob}};
    kwargs...
) where T
    nx, ny = size(first(video))
    d = min(nx, ny)
    nt = length(video)
    p = Dict(
        :maxdist_startvalue => d / 25,
        :maxdist_range => 1:(d / 10),
        :memory_startvalue => 1,
        :memory_range => 0:10,
        :minlife_startvalue => 5,
        :minlife_range => 1:nt,
        :g0_startvalue => 0,
        :g0_range => 0:0.1:2,
        :g2_startvalue => 0,
        :g2_range => 0:0.1:2,
        :w_range => 0:0.05:1,
        :w_startvalue => 0.5,
        :k_range => 1:10,
        :k_startvalue => 1,
        # figure
        :size => (1200, 900),
        :fontsize => 28,
        :fps => 20,
        # blobs
        :blobradius => 4,
        # tracks
        :colormap => :rainbow2,
        :linewidth => 1,
        :maxlength => nt
    )
    # replace parameters if set from kwargs
    for k in keys(kwargs)
        !haskey(p, k) && continue
        p[k] = kwargs[k]
    end
    # initialize
    fig = Figure(; size=p[:size], fontsize=p[:fontsize])
    # controls on bottom
    slider_maxdist = Slider(fig[4,1:4], range=p[:maxdist_range], startvalue=p[:maxdist_startvalue])
    slider_memory = Slider(fig[5,1:4], range=p[:memory_range], startvalue=p[:memory_startvalue])
    slider_minlife = Slider(fig[6,1:4], range=p[:minlife_range], startvalue=p[:minlife_startvalue])
    slider_w = Slider(fig[7,1:4], range=p[:w_range], startvalue=p[:w_startvalue])
    slider_k = Slider(fig[8,1:4], range=p[:k_range], startvalue=p[:k_startvalue])
    label_maxdist = Label(fig[4,5];
        text=@lift("maxdist: "*string($(slider_maxdist.value))),
        halign=:left
    )
    label_memory = Label(fig[5,5];
        text=@lift("memory: "*string($(slider_memory.value))),
        halign=:left
    )
    label_minlife = Label(fig[6,5];
        text=@lift("minlife: "*string($(slider_minlife.value))),
        halign=:left
    )
    label_w = Label(fig[7,5];
        text=@lift("pred weight: "*string($(slider_w.value))),
        halign=:left
    )
    label_w = Label(fig[8,5];
        text=@lift("pred span: "*string($(slider_k.value))),
        halign=:left
    )
    # controls on right
    slider_g0 = Slider(fig[2,4], range=p[:g0_range], startvalue=p[:g0_startvalue])
    slider_g2 = Slider(fig[3,4], range=p[:g2_range], startvalue=p[:g2_startvalue])
    label_g0 = Label(fig[2,5];
        text=@lift("g0: "*string($(slider_g0.value))),
        halign=:left
    )
    label_g2 = Label(fig[3,5];
        text=@lift("g2: "*string($(slider_g2.value))),
        halign=:left
    )
    # video, blobs, trajectories
    t = Observable(1)
    ax1 = Axis(fig[1:3, 1:3];
        aspect=1,
        title="Tracking preview",
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
    )
    img = @lift(video[$t])
    B = @lift(blobs[$t])
    # initialize to the starting values of the sliders
    # then everything will be updated on changes
    trajectories = Observable(
        track_blobs(blobs;
            maxdist=p[:maxdist_startvalue],
            memory=p[:memory_startvalue],
            minlife=p[:minlife_startvalue],
            g0=p[:g0_startvalue],
            g2=p[:g2_startvalue],
            w=p[:w_startvalue],
            k=p[:k_startvalue],
        )
    )
    global heads = @lift(
        [location(track($t)) for track in trajectories[]]
    )
    N = @lift(length($heads))
    labels = @lift([isnantuple($(heads)[i]) ? "" : string(i) for i in eachindex($heads)])
    labelcolors = @lift(resample_cmap(p[:colormap], $N))
    heatmap!(ax1, img; colormap=:bone)
    poly!(ax1, B, p[:blobradius];
        color=:transparent,
        strokecolor=:pink,
        strokewidth=1,
    )
    colormap = @lift(resample_cmap(p[:colormap], length($trajectories)))
    global s = series!(ax1, trajectories, t, p[:maxlength];
        linewidth=p[:linewidth],
        color=colormap,
    )
    global q = text!(ax1, heads; text=labels, color=labelcolors, fontsize=12)
    obs = @lift([
        $(slider_maxdist.value),
        $(slider_memory.value),
        $(slider_minlife.value),
        $(slider_g0.value),
        $(slider_g2.value),
        $(slider_w.value),
        $(slider_k.value),
    ])
    on(obs) do _
        delete!(ax1.scene, s) # clean scene from old tracks
        delete!(ax1.scene, q)
        t[] = 1 # reset time
        trajectories[] = track_blobs(blobs;
            maxdist=slider_maxdist.value[],
            memory=slider_memory.value[],
            minlife=slider_minlife.value[],
            g0=slider_g0.value[],
            g2=slider_g2.value[],
            w=slider_w.value[],
            k=slider_k.value[],
        )
        N[] = length(trajectories[])
        colormap[] = resample_cmap(p[:colormap], N[])
        global s = series!(ax1, trajectories, t, p[:maxlength];
            linewidth=p[:linewidth],
            color=colormap,
        )
        global heads = @lift(
            [location(track($t)) for track in trajectories[]]
        )
        labels[] = [isnantuple(heads[][i]) ? "" : string(i) for i in 1:N[]]
        colormap[] = resample_cmap(p[:colormap], N[])
        global q = text!(ax1, heads; text=labels, color=labelcolors, fontsize=12)
    end
    display(fig)
    while true
        sleep(1/p[:fps])
        if !isopen(fig.scene)
            break
        end
        if t[] < nt
            t[] += 1
        else
            t[] = 1
        end
    end
end
