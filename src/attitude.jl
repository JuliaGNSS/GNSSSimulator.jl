abstract type AbstractAttitude end

struct NoisyStaticAttitude <: AbstractAttitude
    attitude::RotXYZ{Float64}
    base_attitude::RotXYZ{Float64}
    roll_std::Float64
    pitch_std::Float64
    yaw_std::Float64
end

function NoisyStaticAttitude(attitude::RotXYZ, roll_std, pitch_std, yaw_std)
    NoisyStaticAttitude(attitude, attitude, roll_std, pitch_std, yaw_std)
end

struct DynamicAttitude <: AbstractAttitude
    attitudes::Vector{RotXYZ{Float64}}
    time::typeof(1.0s)
    sample_freq::typeof(1.0Hz)
end

@with_kw struct LinearDynamicAttitude <: AbstractAttitude
    attitude::RotXYZ{Float64}
    Δroll::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    Δpitch::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    Δyaw::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
end

struct NoisyDynamicAttitude <: AbstractAttitude
    attitude::RotXYZ{Float64}
    attitudes::Vector{RotXYZ{Float64}}
    time::typeof(1.0s)
    sample_freq::typeof(1.0Hz)
    roll_std::Float64
    pitch_std::Float64
    yaw_std::Float64
end

function NoisyDynamicAttitude(attitudes::Vector{<:RotXYZ}, time, sample_freq, roll_std, pitch_std, yaw_std)
    NoisyDynamicAttitude(first(attitudes), attitudes, time, sample_freq, roll_std, pitch_std, yaw_std)
end

@with_kw struct NoisyLinearDynamicAttitude <: AbstractAttitude
    attitude::RotXYZ{Float64}
    base_attitude::RotXYZ{Float64}
    Δroll::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    Δpitch::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    Δyaw::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    roll_std::Float64
    pitch_std::Float64
    yaw_std::Float64
end

function NoisyLinearDynamicAttitude(attitude::RotXYZ, Δroll, Δpitch, Δyaw, roll_std, pitch_std, yaw_std)
    NoisyLinearDynamicAttitude(attitude, attitude, Δroll, Δpitch, Δyaw, roll_std, pitch_std, yaw_std)
end

function get_attitude(attitude::A) where A <: AbstractAttitude
    attitude.attitude
end

"""
$(SIGNATURES)

Simulates a static noise-free attitude in radian.
"""
function propagate(attitude::RotXYZ, Δt, rng)
    attitude
end

function get_attitude(attitude::RotXYZ)
    attitude
end

"""
$(SIGNATURES)

Simulates a static noisy attitude in radian.
"""
function propagate(attitude::NoisyStaticAttitude, Δt, rng)
    T = attitude.base_attitude
    next_T = RotXYZ(T.theta1 + randn(rng) * attitude.roll_std, T.theta2 + randn(rng) * attitude.pitch_std, T.theta3 + randn(rng) * attitude.yaw_std)
    NoisyStaticAttitude(next_T, T, attitude.roll_std, attitude.pitch_std, attitude.yaw_std)
end

"""
$(SIGNATURES)

Simulates a time-varying noise-free attitude in radian.
If time index exceeds data length, last available value is returned.
"""
function propagate(attitude::DynamicAttitude, Δt, rng)
    next_time = attitude.time + Δt
    DynamicAttitude(attitude.attitudes, next_time, attitude.sample_freq)
end

function get_attitude(attitude::DynamicAttitude)
    get_sampled_value(attitude.time, attitude.sample_freq, attitude.attitudes)
end

"""
$(SIGNATURES)

Simulates a time-varying noisy attitude in radian.
If time index exceeds data length, last available value is returned.
"""
function propagate(attitude::NoisyDynamicAttitude, Δt, rng)
    next_time = attitude.time + Δt
    T = get_sampled_value(attitude.time, attitude.sample_freq, attitude.attitudes)
    next_T = RotXYZ(T.theta1 + randn(rng) * attitude.roll_std, T.theta2 + randn(rng) * attitude.pitch_std, T.theta3 + randn(rng) * attitude.yaw_std)
    NoisyDynamicAttitude(next_T, attitude.attitudes, next_time, attitude.sample_freq, attitude.roll_std, attitude.pitch_std, attitude.yaw_std)
end

"""
$(SIGNATURES)

Simulates a linearly time-varying noise-free attitude in radian.
"""
function propagate(attitude::LinearDynamicAttitude, Δt, rng)
    next_attitude = RotXYZ(attitude.attitude.theta1 + attitude.Δroll * Δt, attitude.attitude.theta2 + attitude.Δpitch * Δt, attitude.attitude.theta3 + attitude.Δyaw * Δt)
    LinearDynamicAttitude(next_attitude, attitude.Δroll, attitude.Δpitch, attitude.Δyaw)
end

"""
$(SIGNATURES)

Simulates a linearly time-varying noisy attitude in radian.
"""
function propagate(attitude::NoisyLinearDynamicAttitude, Δt, rng)
    T = attitude.base_attitude
    next_T = RotXYZ(T.theta1 + attitude.Δroll * Δt, T.theta2 + attitude.Δpitch * Δt, T.theta3 + attitude.Δyaw * Δt)
    noisy_next_T = RotXYZ(next_T.theta1 + randn(rng) * attitude.roll_std, next_T.theta2 + randn(rng) * attitude.pitch_std, next_T.theta3 + randn(rng) * attitude.yaw_std)
    NoisyLinearDynamicAttitude(noisy_next_T, next_T, attitude.Δroll, attitude.Δpitch, attitude.Δyaw, attitude.roll_std, attitude.pitch_std, attitude.yaw_std)
end
