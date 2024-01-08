# API

## [Blobs](@id Blobs)
The `detect_blobs` function processes an image to identify relevant features
in the form of blobs, which are represented by the `Blob` type.
```@docs
AbstractBlob
Blob
```

The properties of a `Blob` can be accessed through the following functions.
```@docs
location
location_raw
scale
amplitude
radius
zeroth_moment
second_moment
intensity_map
```

## [Cost functions](@id Cost)
`blobtracking` links blob into trajectories by finding the optimal solution
to the minimum cost flow problem defined by the keyword `cost` function.
Two cost functions are provided out of the box:
```@docs
QuadraticCost
PCost
```

However, arbitrary cost functions can be defined by the user, and they will
"just work" as long as they respect the API. A cost function must:
- take two and only two positional arguments which are subtypes of `AbstractBlob`
- accept any other relevant argument in the form of a keyword argument
- not modify the state of the blobs
- always return a non-negative real number
Cost functions do not have to explicitly take care of maximum allowed distances,
they should only evaluate the cost of linking the two input blobs.

## [Trajectories](@id Trajectories)
The result of the `blobtracking` is a vector of `Trajectory` objects.
```@docs
Trajectory
```
A `Trajectory` contains a vector of blobs and an associated vector
of time points, which are not necessarily contiguous (in fact, when tracking
with `memory>1`, we allow a trajectory to contain gaps).

All the functions of the Blob API can be also called on `Trajectory` objects,
and, they will return a vector with the value of the requested property
at each step along the trajectory.

A `Trajectory` also supports the following functions
```@docs
length
timestamps
startpoint
endpoint
```
