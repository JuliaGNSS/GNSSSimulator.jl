
const EARTH_RADIUS = 6_360_000
const SPEED_OF_LIGHT = 299_792_458

abstract type AbstractDynamicDOA{T <: Real} end
abstract type AbstractDynamicAttitude{T <: Real} end
abstract type AbstractEmitter end
abstract type AbstractJammer <: AbstractEmitter end

struct GNSSSystem{S<:Function, P<:Function}
    f₀::Float64 # Carrier center frequency
    code_f₀::Float64 # Code center frequency
    code_length::Int
    calc_next_code_phase::P
    gen_sampled_code::S
end

struct EmitterInternalStates
    doppler::Float64
    carrier_phase::Float64
    code_phase::Float64
end

@with_kw struct LinearDynamicDOA{T <: Real} <: AbstractDynamicDOA{T}
    init_DOA::Spherical{T}
    Δazimuth_per_s::T = 0.0
    Δelevation_per_s::T = 0.0
end

@with_kw struct LinearDynamicAttitude{T <: Real} <: AbstractDynamicAttitude{T}
    init_attitude::RotXYZ{T} = RotXYZ(0,0,0)
    Δyaw_per_s::T = 0
    Δpitch_per_s::T = 0
    Δroll_per_s::T = 0
end

"""
$(SIGNATURES)

Returns a `GNSSSystem` struct for L1 parameters.
"""
function gpsl1_system()
    calc_sampled_code, calc_next_code_phase = init_gpsl1_codes()
    GNSSSystem(1_575_420_000., 1_023_000., 1023, calc_next_code_phase, calc_sampled_code)
end

"""
$(SIGNATURES)

Returns a `GNSSSystem` struct for L5 parameters.
"""
function gpsl5_system()
    calc_sampled_code, calc_next_code_phase = init_gpsl5_codes()
    GNSSSystem(1_176_450_000., 10_230_000., 10230, calc_next_code_phase, calc_sampled_code)
end

"""
$(SIGNATURES)

Returns the current DOA based on time `t`. Based on the input `doa` this can be either static if it is of type `Spherical` or
dynamic if it is of type `LinearDynamicDOA`.
"""
function current_enu_doa(t, doa::Spherical)
    doa
end

function current_enu_doa(t, doa::LinearDynamicDOA)
    Spherical(1.0, doa.init_DOA.θ + doa.Δazimuth_per_s * t, doa.init_DOA.ϕ + doa.Δelevation_per_s * t)
end

"""
$(SIGNATURES)

Returns the current attitude based on time `t`. Based on the input `attitude` this can be either static if it is of type `RotXYZ` or
dynamic if it is of type `LinearDynamicAttitude`.
"""
function current_attitude(t, attitude::RotXYZ)
    attitude
end

function current_attitude(t, attitude::LinearDynamicAttitude)
    RotXYZ(attitude.init_attitude.theta1 + attitude.Δroll_per_s * t, attitude.init_attitude.theta2 + attitude.Δpitch_per_s * t, attitude.init_attitude.theta3 + attitude.Δyaw_per_s * t)
end

"""
$(SIGNATURES)

Approximates the distance between user and satellite based on the distance from earth center to satellite `distance_from_earth_center`
and the signal direction of arrival `enu_doa`.
"""
function calc_init_sat_user_distance(distance_from_earth_center, enu_doa)
    init_doa = current_enu_doa(0.0, enu_doa)
    EARTH_RADIUS * cos(init_doa.ϕ + π / 2) + sqrt((EARTH_RADIUS * cos(init_doa.ϕ + π / 2))^2 - EARTH_RADIUS^2 + distance_from_earth_center^2)
end

"""
$(SIGNATURES)

Calculate carrier phase based on the distance between user and satellite `sat_user_distance` and frequency `freq`.
"""
function calc_carrier_phase(sat_user_distance, freq)
    mod2pi(2 * π * freq * sat_user_distance / SPEED_OF_LIGHT)
end

"""
$(SIGNATURES)

Calculates the doppler based on elevation. Assumes satellite velocity only in elevation direction. Currently there is only a positive doppler.
"""
function calc_init_doppler(distance_from_earth_center, enu_doa, velocity, center_freq)
    init_doa = current_enu_doa(0.0, enu_doa)
    center_freq / SPEED_OF_LIGHT * velocity * cos(π / 2 - asin(EARTH_RADIUS * sin(π / 2 + init_doa.ϕ) / distance_from_earth_center))
end

"""
$(SIGNATURES)

Calculate doppler based on the relative velocity `relative_velocity` and the frequency `center_freq`.
"""
function doppler(relative_velocity, center_freq)
    center_freq / SPEED_OF_LIGHT * relative_velocity
end

"""
$(SIGNATURES)

Create a random complex noise signal.
"""
function gen_noise(num_ants, num_samples)
    complex.(randn(num_ants, num_samples), randn(num_ants, num_samples)) / sqrt(2)
end
