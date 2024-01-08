export
    location, pixel, scale, amplitude,
    radius, zeroth_moment, second_moment, intensity_map,
    location_raw

"""
    location(blob)
Return the estimated location of the `blob` including subpixel refinement.
The raw location can be obtained with `location_raw`.
"""
location(blob::AbstractBlob) = blob.location
"""
    scale(blob)
Return the estimated scale of the `blob` along each dimension.
"""
scale(blob::AbstractBlob) = blob.σ
"""
    amplitude(blob)
Return the peak amplitude of the `blob`.
"""
amplitude(blob::AbstractBlob) = blob.amplitude
"""
    radius(blob)
Return the equivalent spherical radius of the `blob`.
"""
radius(blob::AbstractBlob{T,S,N}) where {T,S,N} = sqrt(sum(abs2.(scale(blob))))
"""
    zeroth_moment(blob)
Return the zeroth intensity moment of the `blob` image.
"""
zeroth_moment(blob::AbstractBlob) = blob.m0
"""
    second_moment(blob)
Return the second intensity moment of the `blob` image.
"""
second_moment(blob::AbstractBlob) = blob.m2
"""
    intensity_map(blob)
Return the image of the `blob`.
"""
intensity_map(blob::AbstractBlob) = blob.intensity_map

location_raw(blob::Blob) = blob.location_raw
