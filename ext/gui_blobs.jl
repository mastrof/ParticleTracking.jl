function ParticleTracking.explore_blobs(video::AbstractArray{T,3}; kwargs...) where T
    explore_blobs(eachslice(video; dims=3); kwargs...)
end

function ParticleTracking.explore_blobs(
    video::AbstractVector{<:AbstractMatrix},
    size_range = 1:20, size_startvalue = 1,
    thresh_range = 0.01:0.005:0.15, thresh_startvalue = 0.04,
    ## figure
    size = (1200, 900), fontsize = 28,
    ## blobs
    colormap = :rainbow2,
    markersize = 32,
    alpha = 0.4
)
    ## INITIALIZE
    fig = Figure(; size, fontsize)

    # CONTROLS
    slider_frame = Slider(fig[4,1:4], range=eachindex(video), startvalue=1)
    slider_size = Slider(fig[5,1:4], range=size_range, startvalue=size_startvalue)
    slider_thresh = Slider(fig[6,1:4], range=thresh_range, startvalue=thresh_startvalue)

    label_frame = Label(fig[4,5];
        text=@lift("Frame: "*string($(slider_frame.value))),
        halign=:left
    )
    label_size = Label(fig[5,5];
        text=@lift("Blob size: "*string($(slider_size.value))),
        halign=:left
    )
    label_thresh = Label(fig[6,5];
        text=@lift("Threshold: "*string($(slider_thresh.value))),
        halign=:left
    )

    # MAIN IMAGE WITH BLOBS
    ax1 = Axis(fig[1:3,1:3];
        aspect=1, title="Blob detection preview",
        xticksvisible = false, yticksvisible=false,
        xticklabelsvisible = false, yticklabelsvisible = false,
    )
    img = @lift(video[$(slider_frame.value)])
    blobs = @lift(detect_blobs(
        video[$(slider_frame.value)],
        [$(slider_size.value)];
        rthresh = $(slider_thresh.value)
    ))
    color = @lift(eachindex($blobs)) # use blob ids to assign colors
    ids_each_2 = @lift(eachindex($blobs)[1:2:end])

    heatmap!(ax1, img; colormap=:bone)
    scatter!(ax1, blobs; colormap, color, alpha, markersize)

    # SECONDARY PLOTS WITH BLOB STATISTICS
    ax2 = Axis(fig[1, 4:5]; title="Blob amplitudes",
        xlabel = "id", ylabel = "amplitude",
        ytickformat = "{:.2f}",
        xticks = ids_each_2,
        xminorticksvisible=true, xminorgridvisible=true
    )
    amp = @lift(amplitude.($blobs))
    amprange = @lift(extrema($amp))
    on(amp) do A
        ampmin, ampmax = extrema(A)
        n = length(A)
        xlims!(ax2, 0.5, n+0.5)
        ylims!(ax2, 0.0, 1.15*ampmax)
    end
    hlines!(ax2, slider_thresh.value)
    scatter!(ax2, amp; markersize=20, colormap, color)

    ax3 = Axis(fig[2:3, 4:5]; title="Intensity moments",
        xlabel = rich("m", subscript("0")), ylabel = rich("m", subscript("2")),
        xtickformat = "{:.2f}", ytickformat = "{:.2f}"
    )
    getmoments(blob) = Point2f(zeroth_moment(blob), second_moment(blob))
    m = @lift(getmoments.($blobs))
    on(m) do M
        x1,x2 = extrema(first.(M))
        y1,y2 = extrema(last.(M))
        xlims!(ax3, x1*0.75, x2*1.25)
        ylims!(ax3, y1*0.85, y2*1.15)
    end
    scatter!(ax3, m; markersize=20, colormap, color)

    return fig
end
