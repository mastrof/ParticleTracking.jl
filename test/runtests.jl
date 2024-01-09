using ParticleTracking
using Test

@testset "ParticleTracking.jl" begin
    # draw a gaussian blob centered at point p₀ with scale σ and unit peak intensity
    function gaussian_blob(p, p₀, σ)
        exp(-sum(abs2.(p .- p₀)) / (2*σ^2))
    end
    
    @testset "Blob detection" begin
        n = 256
        for D in 2:2
            @testset "D=$D" begin
                p₀ = Tuple(rand(1:n, D))
                σ = 7
                for noise in [0.1, 0.2, 0.3]
                    # generate noisy background
                    img = rand(ntuple(_ -> n, D)...) .* noise
                    # add blob at p₀
                    for p in CartesianIndices(img)
                        img[p] += gaussian_blob(p.I, p₀, σ)
                        img[p] = min(img[p], 1) # cap at 1
                    end
                    # detect blobs
                    blobs = detect_blobs(img, σ-2:σ+2, rthresh=noise)
                    # if rthresh=noise there should be only one detected blob
                    @test length(blobs) == 1
                    blob = first(blobs)
                    @test location_raw(blob) == CartesianIndex(p₀)
                    # subpixel location should not deviate more than half px from raw
                    @test all(isapprox.(location(blob), p₀; atol=1/2))
                    # scale should match the gaussian standard deviation
                    @test scale(blob) == ntuple(_ -> σ, D)
                end
            end
        end
    end

    @testset "Blob Tracking" begin
        n = 256
        for D in 2:2
            # generate video of moving blob
            p₀ = ntuple(_ -> 25, D)
            v = ntuple(_ -> 10, D) ./ sqrt(D)
            nframes = 10
            real_positions = [p₀ .+ v.*dt for dt in 0:nframes-1]
            noise = 0.1
            σ = 4
            video = [rand(ntuple(_ -> n, D)...).*noise for _ in 1:nframes]
            for i in eachindex(video)
                for p in CartesianIndices(video[i])
                    video[i][p] += gaussian_blob(p.I, real_positions[i], σ)
                    video[i][p] = min(video[i][p], 1) # cap at 1
                end
            end
            blobs = map(frame -> detect_blobs(frame, σ-1:σ+1; rthresh=noise), video)
            @test length(blobs) == nframes
            @test length.(blobs) == ones(Int, nframes)
            trajectories = blobtracking(blobs)
            @test length(trajectories) == 1
            traj = first(trajectories)
            @test length(traj) == nframes-1 # currently last frame is skipped
            # verify that trajectory matches the real positions
            @test all([
                all(isapprox.(location(traj)[i], real_positions[i]; atol=1/2))
                for i in eachindex(blobs)[1:end-1]
            ])
            @test timestamps(traj) == 1:nframes-1

            # test the api
            tracked_blobs = first.(blobs)[1:end-1]
            @test location(traj) == location.(tracked_blobs)
            @test amplitude(traj) == amplitude.(tracked_blobs)
            @test scale(traj) == scale.(tracked_blobs)
            @test radius(traj) == radius.(tracked_blobs)
            @test zeroth_moment(traj) == zeroth_moment.(tracked_blobs)
            @test second_moment(traj) == second_moment.(tracked_blobs)
            @test intensity_map(traj) == intensity_map.(tracked_blobs)
        end
    end
end
