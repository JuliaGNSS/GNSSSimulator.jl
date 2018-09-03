"""
$(SIGNATURES)

Simulates a static noise-free attitude in radian.
"""
function sim_attitude(t, data::StaticAttitudes)
    data.attitude
end

"""
$(SIGNATURES)

Simulates a static noisy attitude in radian.
"""
function sim_attitude(t, data::NoisyStaticAttitudes)
    T = sim_attitude(t, StaticAttitudes(data.attitude))
    RotXYZ(T.theta1 + randn() * data.roll_std, T.theta2 + randn() * data.pitch_std, T.theta3 + randn() * data.yaw_std)
end

"""
$(SIGNATURES)

Simulates a time-varying noise-free attitude in radian.
If time index exceeds data length, last available value is returned
"""
function sim_attitude(t, data::DynamicAttitudes)
    index = floor(Int, t * data.sample_freq) + 1
    index < size(data.attitude, 1) ? data.attitude[index] : data.attitude[end]

end

"""
$(SIGNATURES)

Simulates a time-varying noisy attitude in radian.
If time index exceeds data length, last available value is returned.
"""
function sim_attitude(t, data::NoisyDynamicAttitudes)
    T = sim_attitude(t, DynamicAttitudes(data.attitude, data.sample_freq))
    RotXYZ(T.theta1 + randn() * data.roll_std, T.theta2 + randn() * data.pitch_std, T.theta3 + randn() * data.yaw_std)
end

"""
$(SIGNATURES)

Simulates a linearly time-varying noise-free attitude in radian.
"""
function sim_attitude(t, data::LinearDynamicAttitudes)
    RotXYZ(data.init_attitude.theta1 + data.Δroll * t, data.init_attitude.theta2 + data.Δpitch * t, data.init_attitude.theta3 + data.Δyaw * t)
end

"""
$(SIGNATURES)

Simulates a linearly time-varying noisy attitude in radian.
"""
function sim_attitude(t, data::NoisyLinearDynamicAttitudes)
    T = sim_attitude(t, LinearDynamicAttitudes(data.init_attitude, data.Δroll, data.Δpitch, data.Δyaw))
    RotXYZ(T.theta1 + randn() * data.roll_std, T.theta2 + randn() * data.pitch_std, T.theta3 + randn() * data.yaw_std)
end