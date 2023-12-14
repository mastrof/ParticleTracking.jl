export location, scale, amplitude, radius, zeroth_moment, second_moment, intensity_map

location(blob::AbstractBlob) = blob.location
scale(blob::AbstractBlob) = blob.σ
amplitude(blob::AbstractBlob) = blob.amplitude
radius(blob::AbstractBlob{T,S,N}) where {T,S,N} = scale(blob) .* sqrt(N)
zeroth_moment(blob::AbstractBlob) = blob.m0
second_moment(blob::AbstractBlob) = blob.m2
intensity_map(blob::AbstractBlob) = blob.intensity_map
