# # Tutorial

#=
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
We can manipulate the `size` of the blobs we want to detect and set an
amplitude `threshold` to filter out spurious low-intensity spots which
are wrongly detected as objects.
We can also zoom in and move the image around, as well as check the
quality of the detection on different frames.

On the side, we see a scatter plot showing the zeroth (`m₀`)
and second (`m₂`) intensity moments of the detected
blobs (using a batch of 20 frames to provide some meaningful statistics).
The values of these moments can be used for further filtering or discrimination.
=#

video = Float16.(download("xxx"))
video .-= minimum(video)
video ./= maximum(video)

gui_blobs(video)

#=
After exploring, we can set the desired parameters and actually perform the detection.
=#

blobsize = 2:4
rthresh = 0.05
blobs = detect_blobs(video, blobsize; rthresh)

#=

=#
