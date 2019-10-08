abstract type AbstractSatellite{S <: AbstractGNSSSystem} <: AbstractEmitter end

struct SatellitePhase{T}
    carrier::T
    code::T
end

struct ConstantDopplerSatellite{
    S <: AbstractGNSSSystem,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
} <: AbstractSatellite{S}
    prn::Int
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    code_phase::Float64
    amplitude::Float64
    exists::E
    doa::D
end

function ConstantDopplerSatellite(
        ::Type{S},
        prn::Int;
        carrier_doppler = NaN*Hz,
        carrier_phase = NaN,
        code_phase = NaN,
        amplitude = NaN,
        cn0 = 45dBHz,
        exists::E = true,
        doa::D = SVector(0,0,1),
        velocity = 14_000.0m / 3.6s,
        distance_from_earth_center = 26_560_000.0m,
        n0 = 1/Hz
    ) where {S <: AbstractGNSSSystem,  D <: Union{SVector{3}, AbstractDOA}, E <: Union{Bool, AbstractExistence}}
    if isnan(carrier_doppler)
        carrier_doppler = calc_doppler(distance_from_earth_center, doa, velocity, get_center_frequency(S))
    end
    if isnan(carrier_phase)
        sat_user_distance = calc_sat_user_distance(distance_from_earth_center, doa)
        carrier_phase = calc_carrier_phase(sat_user_distance, get_center_frequency(S) + carrier_doppler)
    end
    if isnan(code_phase)
        sat_user_distance = calc_sat_user_distance(distance_from_earth_center, doa)
        code_doppler = carrier_doppler * get_code_center_frequency_ratio(S)
        code_phase = calc_code_phase(sat_user_distance, get_code_frequency(S) + code_doppler, get_code_length(S))
    end
    if isnan(amplitude)
        amplitude = calc_amplitude_from_cn0(cn0, n0)
    end
    ConstantDopplerSatellite{S, D, E}(prn, carrier_doppler, carrier_phase, code_phase, amplitude, exists, doa)
end

"""
$(SIGNATURES)

Propagates the satellite state over time
"""
function propagate(sat::ConstantDopplerSatellite{S, D, E}, intermediate_frequency, Δt, rng) where {S, D, E, T}
    code_delta = upreferred((get_code_frequency(S) + get_code_doppler(sat)) * Δt)
    code_phase = mod(get_code_phase(sat) + code_delta, get_code_length(S))
    carrier_delta = 2π * upreferred((intermediate_frequency + get_carrier_doppler(sat)) * Δt)
    carrier_phase = mod2pi(get_carrier_phase(sat) + carrier_delta + π) - π
    amplitude = propagate(sat.amplitude, Δt, rng)
    exists = propagate(sat.exists, Δt, rng)
    doa = propagate(sat.doa, Δt, rng)
    ConstantDopplerSatellite{S, D, E}(sat.prn, sat.carrier_doppler, carrier_phase, code_phase, amplitude, exists, doa)
end

function fast_propagate(phase::SatellitePhase, sat::ConstantDopplerSatellite{S}, intermediate_frequency, Δt) where {S <: AbstractGNSSSystem}
    code_length = get_code_length(S)
    code_delta = upreferred((get_code_frequency(S) + get_code_doppler(sat)) * Δt)
    code_phase = phase.code + code_delta
    code_phase = code_phase - (code_phase > code_length) * code_length
    carrier_delta = 2π * upreferred((intermediate_frequency + get_carrier_doppler(sat)) * Δt)
    carrier_phase = phase.carrier + carrier_delta
    carrier_phase = carrier_phase - (carrier_phase > π) * 2π
    SatellitePhase(carrier_phase, code_phase)
end

Base.@propagate_inbounds function get_signal(phase::SatellitePhase, sat::ConstantDopplerSatellite{S}, steer_vec::V, rng) where {S <: AbstractGNSSSystem, T <: Union{Float32, Float64}, V <: Union{SVector{N, Complex{T}}, Complex{T}, T} where N}
    temp = get_existence(sat) *
        T(get_amplitude(sat)) *
        get_code_unsafe(S, phase.code, sat.prn) *
        GNSSSignals.cis_vfast(T(phase.carrier))
    steer_vec * temp
end

@inline get_carrier_doppler(sat::AbstractSatellite) = sat.carrier_doppler
@inline get_code_doppler(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem =
    sat.carrier_doppler * get_code_center_frequency_ratio(S)
@inline get_carrier_phase(sat::AbstractSatellite) = sat.carrier_phase
@inline get_code_phase(sat::AbstractSatellite) = sat.code_phase
@inline get_prn(sat::AbstractSatellite) = sat.prn
@inline get_phase(sat::AbstractSatellite) =
    SatellitePhase(get_carrier_phase(sat), get_code_phase(sat))
@inline get_gnss_system(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem = S

"""
$(SIGNATURES)

Calculate amplitude of a signal based on the carrier-to-noise-density-ratio (CN0) `cn0` in [dB-Hz] and the frequency bandwidth `bandwidth` in [Hz], assumes noise power to be 1.
"""
function calc_amplitude_from_cn0(cn0, n0)
    sqrt(linear(cn0) * n0)
end

"""
$(SIGNATURES)

Approximates the distance between user and satellite based on the distance from earth center to satellite `distance_from_earth_center`
and the signal direction of arrival `enu_doa`.
"""
function calc_sat_user_distance(distance_from_earth_center, doa)
    init_doa = cart2sph(get_doa(doa))
    EARTH_RADIUS * cos(init_doa.ϕ + π / 2) + sqrt((EARTH_RADIUS * cos(init_doa.ϕ + π / 2))^2 - EARTH_RADIUS^2 + distance_from_earth_center^2)
end

"""
$(SIGNATURES)

Calculate carrier phase based on the distance between user and satellite `sat_user_distance` and frequency `freq`.
"""
function calc_carrier_phase(sat_user_distance, freq)
    mod2pi(2π * freq * sat_user_distance / SPEED_OF_LIGHT)
end

"""
$(SIGNATURES)

Calculate code phase based on the distance between satellite and user, the code frequency `freq` and the code length `code_length`.
"""
function calc_code_phase(sat_user_distance, freq, code_length)
    mod(convert(Float64, (freq * sat_user_distance / SPEED_OF_LIGHT)), code_length)
end

"""
$(SIGNATURES)

Calculates the Doppler based on elevation. Assumes satellite velocity only in elevation direction. Currently there is only a positive Doppler.
"""
function calc_doppler(distance_from_earth_center, doa, velocity, center_freq)
    doa_sph = cart2sph(get_doa(doa))
    center_freq / SPEED_OF_LIGHT * velocity * cos(π / 2 - asin(EARTH_RADIUS * sin(π / 2 + doa_sph.ϕ) / distance_from_earth_center))
end
