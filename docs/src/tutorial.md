```@meta
EditURL = "tutorial.jl"
```

# Tutorial

We will use a sample video to test detection and tracking.
The video is available here: [xxx]; we will load it with the
Download.jl package.
We will walk through the entire pipeline, showcasing most of the
package features and how to use them to get the best tracking results.

It is already preprocessed.
For convenience, we want to convert the pixel values to Float16 or Float32
values, and normalize these values in the range [0,1].

````@example tutorial
using Downloads, TiffImages, CairoMakie
fpath = joinpath(pwd(), "sample_video.tif")
Downloads.download(
    "https://github.com/mastrof/ParticleTracking.jl/raw/main/docs/src/sample_video.tif",
    fpath
)
video = TiffImages.load(fpath) .|> Float16
rm(fpath) ## clean up after loading the video
video .-= minimum(video)
video ./= maximum(video)
````

Before performing the full detection, you can explore the data
to find the optimal parameters.

By loading GLMakie, the `explore_blobs` function is made available.
When calling `explore_blobs(video)` you will see an interactive window
where you can move between the different frames of the video and tune
the parameters used for the blob detection: the blob size (i.e. the scale
used in the LoG filtering) and the amplitude threshold.
You can also zoom in and move the image around.

On the side, you will also see two scatter plots, showing the amplitude
of each detected blob, and their zeroth (`m₀`)
and second (`m₂`) intensity moments.


After exploring, we can set the desired parameters and actually perform the detection.
We will limit the blob sizes between 2 to 4 pixels, and set the intensity
threshold to 0.05.

The detection is performed with the `detect_blobs` function, which has to
be applied individually to each frame of `video`.

````@example tutorial
using ParticleTracking

blobsize = 2:4
rthresh = 0.05
blobs = [detect_blobs(img, blobsize; rthresh) for img in eachslice(video; dims=3)]
````

`blobs` is a vector where each element corresponds to a single frame of `video`.
Each element is itself a vector of `Blob`s, an object which contains all the
information about the blob; its position in the image, its size and other properties.
The `Blob` API is described [here](XXX).

A typical step after the initial detection is to look at the zeroth and second
intensity moments of the blobs, which can provide information about differences
between the detected objects.

````@example tutorial
m0 = map.(zeroth_moment, blobs)
m2 = map.(second_moment, blobs)
scatter(vcat(m0...), vcat(m2...); alpha=0.5, axis=(xlabel="m0", ylabel="m2"))
````

The zeroth moment (`m0`) is proportional to the total intensity of the blob;
the tight cluster in the bottom left (`0.5<m0<0.7` and `1.6<m2<1.7`), which
is well separated from the rest of the datapoints,
is then likely to identify spurious detections resulting from pixel shot noise.
If this is the case, it would be good to filter out these blobs
before moving on to the tracking.
A clean set of blobs without spurious detections can drastically
improve the tracking quality, but excessive filtering can have the
opposite effects.

To verify our hypothesis, i.e., low-`m0` low-`m2` blobs are spurious detections,
we have to visually check what these blobs correspond to in the image.
We can then separate the two populations and highlight them in the images
with different colors.

````@example tutorial
maybe_spurious = [
    findall(@. (0.5 < m0[i] < 0.7) && (1.6 < m2[i] < 1.7))
    for i in eachindex(blobs)
]
A = @. view(blobs, maybe_spurious)
B = @. view(blobs, setdiff(eachindex(blobs), maybe_spurious))

# visualize the two populations in frame t
let t = 1
    fig, ax, = heatmap(video[:,:,t], colormap=:bone)
    # blobs can be directly plotted in Makie with scatter
    scatter!(ax, A[t], color=:red, alpha=0.3, markersize=30)
    scatter!(ax, B[t], color=:green, alpha=0.3, markersize=30)
    fig
end
````

The red marker (as can be confirmed by trying different frames `t`), is always
associated to a specific object in the image, but we clearly see it as quite bright!
Why is it different from the others then?
By zooming in, we can see that it is a *single* bright pixel, whereas the other blobs,
while still extremely small, correspond to a larger bright area.

It turns out that the camera I used to collect this video has a dead pixel, which
is always showing up as illuminated.
That is what our moments-based filtering has allowed us to identify.

It is also possible to directly access the `intensity_map` of a blob, i.e. the set of
pixels which have been identified as a unique blob.
We can, for instance, plot the intensity_map of this spurious blob along side that
of another blob (`j`) in the other population.
We see that "real" blobs are generally brighter and/or more extended than the spurious
one.

````@example tutorial
let t=1, I_A=intensity_map(A[t][1]), j=1, I_B=intensity_map(B[t][j])
    M = max(maximum(I_A), maximum(I_B))
    xA, yA = collect.(axes(I_A))
    zA = parent(I_A)
    fig, ax1 = heatmap(xA, yA, zA; colorrange=(0,M), axis=(title="Spurious",))
    ax2 = Axis(fig[1,2], title="Real")
    xB, yB = collect.(axes(I_B))
    zB = parent(I_B)
    heatmap!(ax2, xB, yB, zB; colorrange=(0,M))
    fig
end
````

Now we can safely filter out the spurious detections and move on to the tracking
using the `track_blobs` function.

In principle, we don't need to do anything special and we can just call
`track_blobs(blobs)`, but we will see the results are much less than optimal.

````@example tutorial
# reassign blobs to the subpopulation without spurious detections
blobs = deepcopy(B)
# do the tracking
trajectories = track_blobs(blobs)
# visualize plotting trajectories one by one
let
    fig = Figure()
    ax = Axis(fig[1,1])
    for track in trajectories
        lines!(ax, track)
    end
    fig
end
````

We see that the algorithm tends to connect points that are far away from each other.
Indeed, since the tracking is based on finding the "minimum cost flow" between
successive frames, it always tries to connect points in order to not interrupt
any trajectory. This leads, however, to unreasonable connections.
The first thing to do, therefore, is to define the maximum distance that an object
can cover during a single frame: any connection between points at a distance
larger than this threshold will be less favorable than to just interrupt
the trajectory.

The maximum distance can be set through the `maxdist` keyword.
In our case, after visual inspection it seems reasonable to set
this limit to 60 pixels

````@example tutorial
maxdist = 60
trajectories = track_blobs(blobs; maxdist)
let
    fig = Figure()
    ax = Axis(fig[1,1])
    for track in trajectories
        lines!(ax, track)
    end
    fig
end
````

This is already much better.
But we can do more.

````@example tutorial
scatter(length.(trajectories))
````

A quick look at the length of the trajectories we obtained
shows that we have a few very long trajectories (basically spanning the entire video)
a few medium-length ones and then a bunch of trajectories only 1-frame long.

A 1-frame trajectory is basically a blob that has been detected but not connected to
anything. This can happen for various reasons.
If it's spurious detections that survived our filtering but didn't produce
any connection, that's perfect, we will just filter these trajectories out.
But they may also be real objects which, somehow, were not consistently detected
across successive frames and hence did not produce a trajectory.

In fact, our algorithm has been only looking at immediately successive frames,
but we can also account for the fact that an object might disappear from the
field of view for a few frames, and that we may still be able to recognize it
if it appears again.

Connections across frame gaps are controlled by the `memory` keyword.
By default it is set to 1, meaning that only the immediate next frame is
investigated for a connection. But increasing its value allows us to close
gaps due to objects disappearing from the field of view.
Allowing too much `memory` will, however, have detrimental effects:
the objects are allowed to travel `maxdist` *every frame*, therefore large
`memory` values will allow connections between objects far away in space and time.

````@example tutorial
maxdist = 60
memory = 5
trajectories = track_blobs(blobs; maxdist, memory)
let
    fig = Figure()
    ax = Axis(fig[1,1])
    for track in trajectories
        lines!(ax, track)
    end
    fig
end
````

Compared to the previous memoryless trial, the pink trajectory in the bottom
has been extended; the three nearby blue, green, light-blue segments now constitute
a single, longer trajectory; the two orange and pink trajectories at the top
are now joined into a single one; the brownian trajectories on the left
have not been modified since they were already optimal.

The `track_blobs` function has various additional parameters that can be
tuned to optimize the outcome.

All of them can be explored through an interactive GUI (available via GLMakie)
called by `explore_tracking(video, blobs; kwargs...)`.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

