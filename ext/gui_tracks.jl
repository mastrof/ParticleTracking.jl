isnantuple(t::Tuple) = any(isnan, t)

function in_roi(blob, x1, x2, y1, y2)
    x, y = location(blob)
    (x1 <= x <= x2) && (y1 <= y <= y2)
end

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
        :fontsize => 20,
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
    variables_bottom = ["maxdist", "memory", "minlife", "w", "k"]
    sliders_bottom = Dict(
        variables_bottom .=> map(enumerate(variables_bottom)) do (i, lab)
            var_range = p[Symbol("$(lab)_range")]
            var_start = p[Symbol("$(lab)_startvalue")]
            Slider(fig[3+i,1:4];
                range=var_range, startvalue=var_start,
                # update_while_dragging=false,
            )
        end
    )
    labels_bottom = map(enumerate(variables_bottom)) do (i, lab)
        Label(fig[3+i,5];
            text=@lift("$lab: " * string($(sliders_bottom[lab].value))),
            halign=:left
        )
    end
    # controls on right
    variables_right = ["g0", "g2"]
    sliders_right = Dict(
        variables_right .=> map(enumerate(variables_right)) do (i, lab)
            var_range = p[Symbol("$(lab)_range")]
            var_start = p[Symbol("$(lab)_startvalue")]
            Slider(fig[1+i,4];
                range=var_range, startvalue=var_start,
                # update_while_dragging=false,
            )
        end
    )
    labels_right = Dict(
        variables_right .=> map(enumerate(variables_right)) do (i, lab)
            Label(fig[1+i,5];
                text=@lift("$lab: " * string($(sliders_right[lab].value))),
                halign=:left
            )
        end
    )
    # video, blobs, trajectories
    gl = GridLayout(fig[1:3,1:3])
    t = Observable(1)
    ax1 = Axis(gl[1,1];
        aspect=1,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
    )
    img = @lift(video[$t])
    Nx = size(img[], 2)
    Ny = size(img[], 1)
    xlims!(ax1, 1, Nx)
    ylims!(ax1, 1, Ny)
    # ROI sliders
    roi_x = IntervalSlider(gl[2,1];#fig[3,1:3];
        range=(axes(img[], 2).-1)./Nx,
        startvalues=(0.25, 0.75),
        valign=:bottom,
    )
    roi_y = IntervalSlider(gl[1,2];#fig[1:3,3];
        range=(axes(img[], 1).-1)./Ny,
        startvalues=(0.25, 0.75),
        halign=:right,
        horizontal=false
    )
    colsize!(gl, 1, Relative(0.86))
    colsize!(gl, 2, Relative(0.002))
    rowsize!(gl, 1, Relative(0.998))
    rowsize!(gl, 2, Relative(0.002))
    x1 = @lift($(roi_x.interval)[1]*Nx)
    x2 = @lift($(roi_x.interval)[2]*Nx)
    y1 = @lift($(roi_y.interval)[1]*Ny)
    y2 = @lift($(roi_y.interval)[2]*Ny)
    roi = @lift(Rect(Point2f.([
        ($x1, $y1),
        ($x2, $y2),
    ])))
    plt_roi = poly!(ax1, roi;
        color=:transparent,
        strokewidth=0.1,
        strokecolor=:white,
        linestyle=:dash
    )
    blobs_roi = @lift([
        filter(b -> in_roi(b, $x1, $x2, $y1, $y2), blobs[i])
        for i in eachindex(blobs)
    ])
    B = @lift($(blobs_roi)[$t])
    # initialize to the starting values of the sliders
    # then everything will be updated on changes
    trajectories = Observable(
        track_blobs(blobs_roi[];
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
    plt_img = heatmap!(ax1, img; colormap=:bone)
    plt_blobs = poly!(ax1, B, p[:blobradius];
        color=:transparent,
        strokecolor=:pink,
        strokewidth=1,
    )
    colormap = @lift(resample_cmap(p[:colormap], length($trajectories)))
    global plt_tracks = series!(ax1, trajectories, t, p[:maxlength];
        linewidth=p[:linewidth],
        color=colormap,
    )
    global plt_heads = text!(ax1, heads; text=labels, color=labelcolors, fontsize=12)
    obs = @lift([
        $(sliders_bottom["maxdist"].value),
        $(sliders_bottom["memory"].value),
        $(sliders_bottom["minlife"].value),
        $(sliders_bottom["w"].value),
        $(sliders_bottom["k"].value),
        $(sliders_right["g0"].value),
        $(sliders_right["g2"].value),
    ])
    # on(obs) do
    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.enter
            delete!(ax1.scene, plt_tracks) # clean scene from old tracks
            delete!(ax1.scene, plt_heads)
            t[] = 1 # reset time
            x1[] = roi_x.interval[][1]*Nx
            x2[] = roi_x.interval[][2]Nx
            y1[] = roi_y.interval[][1]*Ny
            y2[] = roi_y.interval[][2]*Ny
            blobs_roi[] = [
                filter(b -> in_roi(b, x1[], x2[], y1[], y2[]), blobs[i])
                for i in eachindex(blobs)
            ]
            trajectories[] = track_blobs(blobs_roi[];
                maxdist=sliders_bottom["maxdist"].value[],
                memory=sliders_bottom["memory"].value[],
                minlife=sliders_bottom["minlife"].value[],
                g0=sliders_right["g0"].value[],
                g2=sliders_right["g2"].value[],
                w=sliders_bottom["w"].value[],
                k=sliders_bottom["k"].value[],
            )
            N[] = length(trajectories[])
            colormap[] = resample_cmap(p[:colormap], N[])
            global plt_tracks = series!(ax1, trajectories, t, p[:maxlength];
                linewidth=p[:linewidth],
                color=colormap,
            )
            global heads = @lift(
                [location(track($t)) for track in trajectories[]]
            )
            labels[] = [isnantuple(heads[][i]) ? "" : string(i) for i in 1:N[]]
            colormap[] = resample_cmap(p[:colormap], N[])
            global plt_heads = text!(ax1, heads; text=labels, color=labelcolors, fontsize=12)
        end
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
