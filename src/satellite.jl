abstract type AbstractSatellite{S <: AbstractGNSSSystem} <: AbstractEmitter end

struct SatellitePhase{T}
    carrier::T
    code::T
end

struct SatellitePhaseWrap
    carrier::Int
    code::Int
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
) where {
    S <: AbstractGNSSSystem,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
}
    if isnan(carrier_doppler)
        carrier_doppler = calc_doppler(
            distance_from_earth_center,
            doa,
            velocity,
            get_center_frequency(S)
        )
    end
    if isnan(carrier_phase)
        sat_user_distance = calc_sat_user_distance(distance_from_earth_center, doa)
        carrier_phase = calc_carrier_phase(
            sat_user_distance,
            get_center_frequency(S) + carrier_doppler
        )
    else
        carrier_phase = carrier_phase / 2π
    end
    if isnan(code_phase)
        sat_user_distance = calc_sat_user_distance(distance_from_earth_center, doa)
        code_doppler = carrier_doppler * get_code_center_frequency_ratio(S)
        code_phase = calc_code_phase(
            sat_user_distance,
            get_code_frequency(S) + code_doppler,
            get_code_length(S) * get_secondary_code_length(S)
        )
    end
    if isnan(amplitude)
        amplitude = calc_amplitude_from_cn0(cn0, n0)
    end
    ConstantDopplerSatellite{S, D, E}(
        prn,
        carrier_doppler,
        carrier_phase,
        code_phase,
        amplitude,
        exists,
        doa
    )
end

"""
$(SIGNATURES)

Propagates the satellite state over time
"""
function propagate(
    sat::ConstantDopplerSatellite{S, D, E},
    intermediate_frequency, Δt, rng
) where {S, D, E, T}
    code_delta = upreferred((get_code_frequency(S) + get_code_doppler(sat)) * Δt)
    code_phase = mod(
        get_code_phase(sat) + code_delta,
        get_code_length(S) * get_secondary_code_length(S)
    )
    carrier_delta = (intermediate_frequency + get_carrier_doppler(sat)) * Δt
    carrier_phase = mod(get_carrier_phase_2pi(sat) + carrier_delta + 0.5, 1) - 0.5
    amplitude = propagate(sat.amplitude, Δt, rng)
    exists = propagate(sat.exists, Δt, rng)
    doa = propagate(sat.doa, Δt, rng)
    ConstantDopplerSatellite{S, D, E}(
        sat.prn,
        sat.carrier_doppler,
        carrier_phase,
        code_phase,
        amplitude,
        exists,
        doa
    )
end

function calc_phase(
    sat::ConstantDopplerSatellite{S},
    index,
    intermediate_frequency,
    sample_frequency
) where {S <: AbstractGNSSSystem}
    carrier_phase = (intermediate_frequency + get_carrier_doppler(sat)) * index /
        sample_frequency + get_carrier_phase_2pi(sat)
    code_phase = (get_code_frequency(S) + get_code_doppler(sat)) * index /
        sample_frequency + get_code_phase(sat)
    SatellitePhase(carrier_phase, code_phase)
end

function init_phase_wrap(sat::ConstantDopplerSatellite)
    SatellitePhaseWrap(0, 0)
end

function update_phase_wrap(
    phase_wrap::SatellitePhaseWrap,
    phase::SatellitePhase,
    sat::ConstantDopplerSatellite{S},
) where {S <: AbstractGNSSSystem}
    code_length = get_code_length(S) * get_secondary_code_length(S)
    code_wrap = phase_wrap.code + ((phase.code - phase_wrap.code) >= code_length) *
        code_length
    carrier_wrap = phase_wrap.carrier + ((phase.carrier - phase_wrap.carrier) >= 0.5) * 1
    SatellitePhaseWrap(carrier_wrap, code_wrap)
end

Base.@propagate_inbounds function get_signal(
    phase::SatellitePhase,
    phase_wrap::SatellitePhaseWrap,
    sat::ConstantDopplerSatellite{S},
    steer_vec::V,
    rng
) where {
    S <: AbstractGNSSSystem,
    T <: Union{Float32, Float64},
    V <: Union{SVector{N, Complex{T}}, Complex{T}, T} where N
}
    temp = get_existence(sat) *
        T(get_amplitude(sat)) *
        get_code_unsafe(S, phase.code - phase_wrap.code, sat.prn) *
        GNSSSignals.cis_vfast(T(2π * phase.carrier - phase_wrap.carrier))
    steer_vec * temp
end

@inline get_carrier_doppler(sat::AbstractSatellite) = sat.carrier_doppler
@inline get_code_doppler(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem =
    sat.carrier_doppler * get_code_center_frequency_ratio(S)
@inline get_carrier_phase(sat::AbstractSatellite) = 2π * sat.carrier_phase
@inline get_carrier_phase_2pi(sat::AbstractSatellite) = sat.carrier_phase
@inline get_code_phase(sat::AbstractSatellite) = sat.code_phase
@inline get_prn(sat::AbstractSatellite) = sat.prn
@inline get_gnss_system(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem = S

"""
$(SIGNATURES)

Calculate amplitude of a signal based on the carrier-to-noise-density-ratio (CN0) `cn0` in
[dB-Hz] and the frequency bandwidth `bandwidth` in [Hz], assumes noise power to be 1.
"""
function calc_amplitude_from_cn0(cn0, n0)
    sqrt(linear(cn0) * n0)
end

"""
$(SIGNATURES)

Approximates the distance between user and satellite based on the distance from earth center
to satellite `distance_from_earth_center` and the signal direction of arrival `enu_doa`.
"""
function calc_sat_user_distance(distance_from_earth_center, doa)
    init_doa = cart2sph(get_doa(doa))
    EARTH_RADIUS * cos(init_doa.ϕ + π / 2) +
        sqrt((EARTH_RADIUS * cos(init_doa.ϕ + π / 2))^2 - EARTH_RADIUS^2 +
        distance_from_earth_center^2)
end

"""
$(SIGNATURES)

Calculate carrier phase based on the distance between user and satellite `sat_user_distance`
and frequency `freq`.
"""
function calc_carrier_phase(sat_user_distance, freq)
    mod(freq * sat_user_distance / SPEED_OF_LIGHT + 1, 2) - 1
end

"""
$(SIGNATURES)

Calculate code phase based on the distance between satellite and user, the code frequency
`freq` and the code length `code_length`.
"""
function calc_code_phase(sat_user_distance, freq, code_length)
    mod(convert(Float64, (freq * sat_user_distance / SPEED_OF_LIGHT)), code_length)
end

"""
$(SIGNATURES)

Calculates the Doppler based on elevation. Assumes satellite velocity only in elevation
direction. Currently there is only a positive Doppler.
"""
function calc_doppler(distance_from_earth_center, doa, velocity, center_freq)
    doa_sph = cart2sph(get_doa(doa))
    center_freq / SPEED_OF_LIGHT * velocity * cos(π / 2 - asin(EARTH_RADIUS *
        sin(π / 2 + doa_sph.ϕ) / distance_from_earth_center))
end
