abstract type AbstractSatellite{S <: AbstractGNSSSystem} <: AbstractEmitter end

struct ConstantDopplerSatellite{
    S <: AbstractGNSSSystem,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
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
    amplitude = calc_amplitude_from_cn0(cn0, n0)
    ConstantDopplerSatellite{S, D, E}(prn, carrier_doppler, carrier_phase, code_phase, amplitude, exists, doa)
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

"""
$(SIGNATURES)

Propagates the satellite state over time
"""
function propagate(sat::ConstantDopplerSatellite{S, D, E}, Δt) where {S <: AbstractGNSSSystem,  D <: Union{SVector{3}, AbstractDOA}, E <: Union{Bool, AbstractExistence}}
    carrier_phase = 2π * get_carrier_doppler(sat) * Δt + get_carrier_phase(sat)
    carrier_phase -= (carrier_phase > 2π) * 2π
    code_phase = (get_code_frequency(S) + get_code_doppler(sat)) * Δt + get_code_phase(sat)
    code_length = get_code_length(S)
    code_phase -= (code_phase >= code_length) * code_length
    amplitude = propagate(sat.amplitude, Δt)
    exists = propagate(sat.exists, Δt)
    doa = propagate(sat.doa, Δt)
    ConstantDopplerSatellite{S, D, E}(sat.prn, sat.carrier_doppler, carrier_phase, code_phase, amplitude, exists, doa)
end

function get_signal(sat::AbstractSatellite{S}, attitude, manifold) where S <: AbstractGNSSSystem
    get_existence(sat) * (get_steer_vec(manifold, get_doa(sat), attitude) * (
        cis(get_carrier_phase(sat)) *
        get_code_unsafe(S, get_code_phase(sat), get_prn(sat)) *
        get_amplitude(sat)
    ))
end

get_carrier_doppler(sat::AbstractSatellite) = sat.carrier_doppler
get_code_doppler(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem = sat.carrier_doppler * get_code_center_frequency_ratio(S)
get_carrier_phase(sat::AbstractSatellite) = sat.carrier_phase
get_code_phase(sat::AbstractSatellite) = sat.code_phase
get_prn(sat::AbstractSatellite) = sat.prn

"""
$(SIGNATURES)

Calculate amplitude of a signal based on the carrier-to-noise-density-ratio (CN0) `cn0` in [dB-Hz] and the frequency bandwidth `bandwidth` in [Hz], assumes noise power to be 1.
"""
function calc_amplitude_from_cn0(cn0, n0)
    sqrt(linear(cn0) * n0)
end
