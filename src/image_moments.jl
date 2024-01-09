"""
    image_moments(location, σ, img)
Evaluate the zeroth and second intensity moment of a blob in `img`
identified as a a region of size `σ` around `location`.
Returns a tuple with a view of the blob and the two moments.
"""
function image_moments(location::CartesianIndex{N}, σ, img::AbstractArray{T,N}) where {T,N}
    x = location.I
    σ2 = σ.^2
    m0 = zero(float(T))
    m2 = zero(float(T))
    # define rectangle around the blob
    view_range = ntuple(i -> (x[i]-σ[i]):(x[i]+σ[i]), N)
    img_padded = PaddedView(0.0, img, view_range)
    I = zero(img_padded)
    for c in CartesianIndices(img_padded)
        delta2 = (c.I .- x).^2
        if sum(delta2./σ2) <= 1 # N-sphere of radius σ
            A = img_padded[c]
            I[c] = A
            m0 += A
            m2 += sum(delta2)*A
        end
    end
    m2 /= m0
    return I, m0, m2
end
