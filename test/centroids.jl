using ParticleTracking
using Test

@testset "Centroids" begin
    nx, ny = 10, 10
    img = zeros(nx, ny)
    img[2, 2] = 1
    img[6:7, 4:6] .= 1
    cc = connected_components(img)
    centroids = component_centroids(cc) # image coordinates
    @test centroids == [(2.0, 2.0), (6.5, 5.0)]
    xs = range(0, 1; length=nx)
    ys = range(1, 2; length=ny)
    centroids = component_centroids(cc, (xs,ys)) # real-space coordinates
    @test (
        all(centroids[1] .≈ (xs[2], ys[2]))
        &&
        all(centroids[2] .≈ ((xs[6]+xs[7])/2, ys[5]))
    )

    nx, ny, nz = 10, 10, 10
    img1 = zeros(nx, ny, nz); img1[2,2,6] = 1
    img2 = zeros(nx, ny, nz); img2[2,3,5:6] = 1
    stack = [img1, img2]
    cc = connected_components(stack)
    centroids = component_centroids(cc) # image coordinates
    @test centroids == [[(2.0, 2.0, 6.0)], [(2.0, 3.0, 5.5)]]
    xs = range(0, 1; length=nx)
    ys = range(0, 2; length=ny)
    zs = range(-0.5, 0.5; length=nz)
    centroids = component_centroids(cc, (xs,ys,zs))
    @test centroids == [[(xs[2], ys[2], zs[6])], [(xs[2], ys[3], (zs[5]+zs[6])/2)]]
end