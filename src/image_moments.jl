function image_moments(location::CartesianIndex{2}, σ, img::AbstractArray{T,2}) where T
    x, y = location.I
    σx, σy = σ
    σx2 = σx^2
    σy2 = σy^2
    m0 = zero(float(T))
    m2 = zero(float(T))
    img_padded = PaddedView(0.0, img, (x-σx:x+σx, y-σy:y+σy))
    I = zero(img_padded)
    for j in -σy:σy, i in -σx:σx
        i2 = i^2
        j2 = j^2
        if i2/σx2 + j2/σy2 <= 1
            A = img_padded[x+i, y+j]
            I[x+i, y+j] = A
            m0 += A
            m2 += (i2 + j2) * A
        end
    end
    #m2 /= (m0 * σx * σy)
    #m0 /= (σx * σy)
    m2 /= m0
    return I, m0, m2
end
