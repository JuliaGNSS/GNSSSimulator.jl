
const EARTH_RADIUS = 6_360_000
const SPEED_OF_LIGHT = 299_792_458

abstract type AbstractDynamicDOA{T <: Real} end
abstract type AbstractDynamicAttitude{T <: Real} end
abstract type AbstractEmitter end
abstract type AbstractJammer <: AbstractEmitter end

struct GNSSSystem
    f₀::Float64 # Carrier center frequency
    code_f₀::Float64 # Code center frequency
    code_length::Int
    calc_next_code_phase::Function
    gen_sampled_code::Function
end

function gpsl1_system()
    calc_sampled_code, calc_next_code_phase = init_gpsl1_codes()
    GNSSSystem(1_575_420_000, 1_023_000, 1023, calc_next_code_phase, calc_sampled_code)
end

function gpsl5_system()
    calc_sampled_code, calc_next_code_phase = init_gpsl5_codes()
    GNSSSystem(1_176_450_000, 10_230_000, 10230, calc_next_code_phase, calc_sampled_code)
end

@with_kw struct LinearDynamicDOA{T <: Real} <: AbstractDynamicDOA{T}
    init_DOA::Spherical{T}
    Δazimuth_per_s::T = 0
    Δelevation_per_s::T = 0
end

@with_kw struct Satellite{T <: Union{Spherical{R}, AbstractDynamicDOA{R}} where R<:Real} <: AbstractEmitter
    svid::Int
    enu_doa::T
    velocity::Float64 = 14_000 / 3.6 # Assume sat velocity only in elevation direction
    CN0::Float64 = 45
    distance_from_earth_center::Float64 = 26_560_000
end

struct CWJammer{T <: Union{Spherical{R}, AbstractDynamicDOA{R}} where R<:Real} <: AbstractJammer
    id::Int
    enu_doa::T
    relative_velocity::Float64
    JNR::Float64
end

@with_kw struct LinearDynamicAttitude{T <: Real} <: AbstractDynamicAttitude{T}
    init_attitude::RotXYZ{T} = RotXYZ(0,0,0)
    Δyaw_per_s::T = 0
    Δpitch_per_s::T = 0
    Δroll_per_s::T = 0
end

function current_enu_doa(t, doa::Spherical)
    doa
end

function current_enu_doa(t, doa::LinearDynamicDOA)
    Spherical(1.0, doa.init_DOA.θ + doa.Δazimuth_per_s * t, doa.init_DOA.ϕ + doa.Δelevation_per_s * t)
end

function current_attitude(t, attitude::RotXYZ)
    attitude
end

function current_attitude(t, attitude::LinearDynamicAttitude)
    RotXYZ(attitude.init_attitude.theta1 + attitude.Δroll_per_s, attitude.init_attitude.theta2 + attitude.Δpitch_per_s, attitude.init_attitude.theta3 + attitude.Δyaw_per_s)
end

function init_sat_distance(distance_from_earth_center, enu_doa)
    doa = current_enu_doa(0.0, enu_doa)
    EARTH_RADIUS * cos(doa.ϕ + π / 2) + sqrt((EARTH_RADIUS * cos(doa.ϕ + π / 2))^2 - EARTH_RADIUS^2 + distance_from_earth_center^2)
end

function carrier_phase(sat_user_distance, freq)
    mod2pi(2 * π * freq * sat_user_distance / SPEED_OF_LIGHT)
end

function code_phase(sat_user_distance, freq, code_length)
    mod(freq * sat_distance / SPEED_OF_LIGHT, code_length)
end

function doppler(distance_from_earth_center, enu_doa, velocity, center_freq) # Assume sat velocity only in elevation direction
    center_freq / SPEED_OF_LIGHT * velocity * cos(π / 2 - asin(EARTH_RADIUS * sin(π / 2 + enu_doa.init_DOA.ϕ) / distance_from_earth_center))
end

function doppler(relative_velocity, center_freq)
    center_freq / SPEED_OF_LIGHT * relative_velocity
end

function calc_amplitude_from_cn0(cn0_dB, bandwidth) # Assumes noise power of 1
    sqrt(10^(cn0_dB / 10) / bandwidth)
end

function calc_amplitude_from_jnr(jnr_dB) # Assumes noise power of 1
    10^(jnr_dB / 20)
end

function gen_noise(num_ants, num_sampels)
    complex.(randn(num_ants, num_samples), randn(num_ants, num_samples)) / sqrt(2)
end
