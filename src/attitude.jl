"""
$(SIGNATURES)

Simulates a static noise-free attitude in radian.
"""
function sim_attitude(t, data::RotXYZ{Float64})
    data
end

"""
$(SIGNATURES)

Simulates a static noisy attitude in radian.
"""
function sim_attitude(t, data::NoisyStaticAttitude)
    T = sim_attitude(t, data.attitude)
    RotXYZ(T.theta1 + randn() * data.roll_std, T.theta2 + randn() * data.pitch_std, T.theta3 + randn() * data.yaw_std)
end

"""
$(SIGNATURES)

Simulates a time-varying noise-free attitude in radian.
If time index exceeds data length, last available value is returned.
"""
function sim_attitude(t, data::DynamicAttitude)
    index = floor(Int, t * data.sample_freq) + 1
    index < size(data.attitude, 1) ? data.attitude[index] : data.attitude[end]

end

"""
$(SIGNATURES)

Simulates a time-varying noisy attitude in radian.
If time index exceeds data length, last available value is returned.
"""
function sim_attitude(t, data::NoisyDynamicAttitude)
    T = sim_attitude(t, DynamicAttitude(data.attitude, data.sample_freq))
    RotXYZ(T.theta1 + randn() * data.roll_std, T.theta2 + randn() * data.pitch_std, T.theta3 + randn() * data.yaw_std)
end

"""
$(SIGNATURES)

Simulates a linearly time-varying noise-free attitude in radian.
"""
function sim_attitude(t, data::LinearDynamicAttitude)
    RotXYZ(data.init_attitude.theta1 + data.Δroll * t, data.init_attitude.theta2 + data.Δpitch * t, data.init_attitude.theta3 + data.Δyaw * t)
end

"""
$(SIGNATURES)

Simulates a linearly time-varying noisy attitude in radian.
"""
function sim_attitude(t, data::NoisyLinearDynamicAttitude)
    T = sim_attitude(t, LinearDynamicAttitude(data.init_attitude, data.Δroll, data.Δpitch, data.Δyaw))
    RotXYZ(T.theta1 + randn() * data.roll_std, T.theta2 + randn() * data.pitch_std, T.theta3 + randn() * data.yaw_std)
end