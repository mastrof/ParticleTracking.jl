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
        :memory_range => 1:10,
        :minlife_startvalue => 5,
        :minlife_range => 1:nt,
        :g0_startvalue => 0,
        :g0_range => 0:0.1:2,
        :g2_startvalue => 0,
        :g2_range => 0:0.1:2,
        # figure
        :size => (1200, 900),
        :fontsize => 28,
        :fps => 20,
        # blobs
        :blobradius => 32,
        :alpha => 0.2,
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
    trajectories = @lift(
        blobtracking(blobs;
            maxdist=$(slider_maxdist.value),
            memory=$(slider_memory.value),
            minlife=$(slider_minlife.value),
            g0=$(slider_g0.value),
            g2=$(slider_g2.value),
        )
    )
    heads = @lift([location(track($t)) for track in $trajectories])
    N = @lift(length($trajectories))
    labels = @lift(string.(1:$N))
    labelcolors = @lift(resample_cmap(p[:colormap], $N))
    heatmap!(ax1, img; colormap=:bone)
    scatter!(ax1, B;
        color=:pink,
        alpha=p[:alpha],
        markersize=p[:markersize]
    )
    on(trajectories) do tracks
        notify(heads)
        notify(N)
        notify(labels)
        colormap = resample_cmap(p[:colormap], length(tracks))
        for (i, s) in enumerate(tracks)
            lines!(ax1, s, t, p[:maxlength];
                linewidth=p[:linewidth],
                color=colormap[i]
            )
        end
        text!(ax1, heads; text=labels, color=labelcolors, fontsize=12)
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


# ## GUI
# fig = Figure(; resolution=(1200,900), fontsize=28)
# stack = proc[1:200]
# allblobs = detect_blobs(stack, 2:4; rthresh=0.05)
#
# # CONTROLS
# slider_frame = Slider(fig[4,1:4], range=eachindex(stack), startvalue=1)
# slider_memory = Slider(fig[5,1:4], range=1:10, startvalue=1)
# slider_maxdist = Slider(fig[6,1:4], range=5:5:75, startvalue=30)
#
# label_frame = Label(fig[4,5];
#     text=@lift("Frame: "*string($(slider_frame.value))),
#     halign=:left
# )
# label_memory = Label(fig[5,5];
#     text=@lift("Memory: "*string($(slider_memory.value))),
#     halign=:left
# )
# label_maxdist = Label(fig[6,5];
#     text=@lift("Max distance: "*string($(slider_maxdist.value))),
#     halign=:left
# )
#
# # MAIN IMAGE WITH BLOBS AND TRACKS
# ax1 = Axis(fig[1:3,1:3];
#     aspect=1, title="Blob detection preview",
#     xticksvisible = false, yticksvisible=false,
#     xticklabelsvisible = false, yticklabelsvisible = false,
# )
# img = @lift(Gray.(stack[$(slider_frame.value)]))
# blobs = @lift(allblobs[$(slider_frame.value)])
# tracks = @lift(blobtracking(allblobs;
#     memory=$(slider_memory.value),
#     maxdist=$(slider_maxdist.value),
#     minlife=5
# ))
#
# heatmap!(ax1, img)
# scatter!(ax1, blobs; markersize=32, alpha=0.2, color=:red)
# on(tracks) do _tracks
#     for s in _tracks
#         lines!(ax1, s, slider_frame.value, 50; linewidth=2)
#     end
# end
#
# # SECONDARY PLOTS WITH BLOB STATISTICS
# ax2 = Axis(fig[1, 4:5]; title="Blob amplitudes",
#     xlabel = "id", ylabel = "amplitude"
# )
# amp = @lift(amplitude.($blobs))
# amprange = @lift(extrema($amp))
# on(amp) do A
#     ampmin, ampmax = extrema(A)
#     n = length(A)
#     xlims!(ax2, 0.5, n+0.5)
#     ylims!(ax2, 0.0, 1.15*ampmax)
# end
# scatter!(ax2, amp)
#
# ax3 = Axis(fig[2:3, 4:5]; title="Intensity moments",
#     xlabel = rich("m", subscript("0")), ylabel = rich("m", subscript("2"))
# )
# getmoments(blob) = Point2f(zeroth_moment(blob), second_moment(blob))
# m = @lift(getmoments.($blobs))
# on(m) do M
#     x1,x2 = extrema(first.(M))
#     y1,y2 = extrema(last.(M))
#     xlims!(ax3, x1*0.75, x2*1.25)
#     ylims!(ax3, y1*0.75, y2*1.25)
# end
# scatter!(ax3, m)
#
#
# fig
