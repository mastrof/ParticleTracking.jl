## GUI
fig = Figure(; resolution=(1200,900), fontsize=28)

# CONTROLS
slider_frame = Slider(fig[4,1:4], range=eachindex(proc), startvalue=1)
slider_size = Slider(fig[5,1:4], range=1:20, startvalue=1)
slider_thresh = Slider(fig[6,1:4], range=0.01:0.01:0.4, startvalue=0.03)

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
img = @lift(Gray.(proc[$(slider_frame.value)]))
blobs = @lift(detect_blobs(
    proc[$(slider_frame.value)],
    [$(slider_size.value)];
    rthresh = $(slider_thresh.value)
))

heatmap!(ax1, img)
scatter!(ax1, blobs; markersize=32, alpha=0.2, color=:red)

# SECONDARY PLOTS WITH BLOB STATISTICS
ax2 = Axis(fig[1, 4:5]; title="Blob amplitudes",
    xlabel = "id", ylabel = "amplitude"
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
scatter!(ax2, amp)

ax3 = Axis(fig[2:3, 4:5]; title="Intensity moments",
    xlabel = rich("m", subscript("0")), ylabel = rich("m", subscript("2"))
)
getmoments(blob) = Point2f(zeroth_moment(blob), second_moment(blob))
m = @lift(getmoments.($blobs))
on(m) do M
    x1,x2 = extrema(first.(M))
    y1,y2 = extrema(last.(M))
    xlims!(ax3, x1*0.75, x2*1.25)
    ylims!(ax3, y1*0.75, y2*1.25)
end
scatter!(ax3, m)


fig
