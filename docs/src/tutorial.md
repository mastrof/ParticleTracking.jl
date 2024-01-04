```@meta
EditURL = "tutorial.jl"
```

# Tutorial

We will use a sample video to test detection and tracking.
The video is available here: [xxx]; we will load it with the
Download.jl package.

It is already preprocessed.
For convenience, we want to convert the pixel values to Float16 or Float32
values, and normalize these values in the range [0,1].
Before performing the full detection, we can explore the data
to find the optimal parameters.

An interactive GUI is provided through GLMakie with the `gui_blobs` function.
In the large main window we see one frame of the video, with the detected blobs
surrounded by red circles.
We can manipulate the `size` of the blobs we want to detect and an
amplitude `threshold` to filter out spurious low-intensity spots which
are wrongly detected as objects.

On the side, we see a scatter plot showing the zeroth (`m₀`)
and second (`m₂`) intensity moments of the detected
blobs (using a batch of 20 frames to provide some meaningful statistics).
The values of these moments can be used for further filtering or discrimination.

````@example tutorial
video = ".."
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

