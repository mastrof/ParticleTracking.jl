export timestamps, startpoint, endpoint
export location, amplitude, scale, radius, zeroth_moment, second_moment, intensity_map
export displacement_net, displacement_gross


"""
    timestamps(trajectory)
Return the times along the `trajectory` at which blobs have been measured.
"""
timestamps(traj::Trajectory) = traj.times
"""
    length(trajectory)
Return the length of the `trajectory`.
"""
Base.length(traj::Trajectory) = length(traj.times)
"""
    startpoint(trajectory)
Return the time at which the `trajectory` started.
"""
startpoint(traj::Trajectory) = location(first(traj.blobs))
"""
    endpoint(trajectory)
Return the time at which the `trajectory` ended.
"""
endpoint(traj::Trajectory) = location(last(traj.blobs))

for f in (location, amplitude, scale, zeroth_moment, second_moment, intensity_map)
    m = parentmodule(f)
    fs = nameof(f)
    @eval function $m.$fs(traj::Trajectory, args...)
        t1, t2 = extrema(timestamps(traj))
        #return OffsetArray($fs.(traj.blobs), t1:t2)
        return $fs.(traj.blobs)
    end
end

@inline function displacement_net(traj::Trajectory)
    x = startpoint(traj)
    y = endpoint(traj)
    sqrt(sum(abs2.(x .- y)))
end

@inline function displacement_gross(traj::Trajectory)
    positions = location(traj)
    d2 = [sum(abs2.(positions[i] .- positions[i.-1])) for i in eachindex(timestamps(traj))[2:end]]
    sum(sqrt.(d2))
end
