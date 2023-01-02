using ParticleTracking
using Test

@testset "Connected components" begin
    @testset "Two dimensions" begin
        nx, ny = 20, 20
        img = zeros(nx, ny)
        # create two connected components
        img[9:11, 9:11] .= 1
        # close to the other one but has different value → not connected
        img[8,9] = 0.5
        cc = connected_components(img)
        c₁ = vec([CartesianIndex(i,j) for (i,j) in Iterators.product(9:11, 9:11)])
        c₂ = [CartesianIndex(8,9)]
        @test length(cc) == 2 && Set(cc) == Set((c₁, c₂))

        nx, ny = 10, 10
        img = zeros(nx, ny)
        # create a connected component
        img[4:6, 4:6] .= 1
        # connectivity only along x direction
        cc = connected_components(img, dims=(1,))
        # there should be three connected components (4:6,4), (4:6,5), (4:6,6)
        real_cc = [CartesianIndex.(4:6,j) for j in 4:6]
        @test length(cc) == 3 && Set(cc) == Set(real_cc)
        # only connect along diagonals with connectivity matrix
        connectivity = BitMatrix([1 0 1; 0 1 0; 1 0 1])
        cc = connected_components(img, connectivity)
        # now there should be two distinct connected components
        c₁ = CartesianIndex.([ (4,4), (6,4), (5,5), (4,6), (6,6) ])
        c₂ = CartesianIndex.([ (5,4), (4,5), (6,5), (5,6) ])
        @test length(cc) == 2 && Set(cc) == Set((c₁, c₂))

        nx, ny = 10, 10
        img = zeros(nx, ny)
        img[5, 5] = 1
        # if bkg=1, the whole image is a connected component except (5,5)
        cc = connected_components(img; bkg=1)
        c₁ = [CartesianIndex(i,j) for i in 1:nx, j in 1:ny if (i,j)≠(5,5)]
        @test length(cc) == 1 && first(cc) == c₁

        # check that dispatch on stacks works as expected
        nx, ny = 10, 10
        img1 = zeros(nx,ny); img1[5,5] = 1
        img2 = zeros(nx,ny); img2[5,6] = 1
        stack = [img1, img2]
        cc = connected_components(stack)
        @test cc[1] == [[CartesianIndex(5,5)]] && cc[2] == [[CartesianIndex(5,6)]]

        img1 = zeros(nx,ny).+0.5; img1[5,5] = 1
        img2 = zeros(nx,ny).+0.5; img2[5,6] = 1
        stack = [img1, img2]
        cc = connected_components(stack; bkg=0.5)
        @test cc[1] == [[CartesianIndex(5,5)]] && cc[2] == [[CartesianIndex(5,6)]]

        img1 = zeros(nx,ny); img1[5,5:6] .= 1
        img2 = zeros(nx,ny); img2[5:6,5] .= 1
        stack = [img1, img2]
        cc = connected_components(stack, dims=(2,)) # only y connections
        @test (
            cc[1] == [CartesianIndex.(5,5:6)]
            &&
            cc[2] == [[CartesianIndex(5,5)], [CartesianIndex(6,5)]]
        )
    end

    # only a few extra tests for safety; if 2D works, 3D should work as well
    @testset "Three dimensions" begin
        nx, ny, nz = 20, 20, 5
        img = zeros(nx, ny, nz)
        # create two connected components
        img[9:11, 9:11, 2:3] .= 1
        # close to the other one but has different value → not connected
        img[8,9,3] = 0.5
        cc = connected_components(img)
        c₁ = vec([CartesianIndex(i,j,k) for (i,j,k) in Iterators.product(9:11,9:11,2:3)])
        c₂ = [CartesianIndex(8,9,3)]
        @test length(cc) == 2 && Set(cc) == Set((c₁, c₂))

        # don't connect along z
        nx, ny, nz = 20, 20, 5
        img = zeros(nx, ny, nz)
        img[9:11, 9:11, 2:3] .= 1
        cc = connected_components(img, dims=(1,2))
        c₁ = vec([CartesianIndex(i,j,2) for (i,j) in Iterators.product(9:11,9:11)])
        c₂ = vec([CartesianIndex(i,j,3) for (i,j) in Iterators.product(9:11,9:11)])
        @test length(cc) == 2 && Set(cc) == Set((c₁, c₂))
    end
end