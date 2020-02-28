abstract type AbstractSatellite{S <: AbstractGNSSSystem, T} <: AbstractEmitter{T} end

struct ConstantDopplerSatellite{
    S <: AbstractGNSSSystem,
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Unitful.Level{Unitful.Decibel, 1Hz, <:Unitful.Quantity{T}}
} <: AbstractSatellite{S, T}
    prn::Int
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    code_phase::Float64
    cn0::C
    exists::E
    doa::D
    carrier_code::StructArray{Complex{Int16},1,NamedTuple{(:re, :im),Tuple{Array{Int16,1},Array{Int16,1}}},Int}
    signal::StructArray{Complex{T},1,NamedTuple{(:re, :im),Tuple{Array{T,1},Array{T,1}}},Int}
end

function ConstantDopplerSatellite(
    ::Type{S},
    prn::Int;
    carrier_doppler = NaN*Hz,
    carrier_phase = NaN,
    code_phase = NaN,
    cn0::C = 45dBHz,
    exists::E = true,
    doa::D = SVector(0,0,1),
    velocity = 14_000.0m / 3.6s,
    distance_from_earth_center = 26_560_000.0m
) where {
    S <: AbstractGNSSSystem,
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Unitful.Level{Unitful.Decibel, 1Hz, <:Unitful.Quantity{T}},
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
    carrier_code = StructArray{Complex{Int16}}(undef, 0)
    signal = StructArray{Complex{T}}(undef, 0)
    ConstantDopplerSatellite{S, T, D, E, C}(
        prn,
        carrier_doppler,
        carrier_phase,
        code_phase,
        cn0,
        exists,
        doa,
        carrier_code,
        signal
    )
end

function gen_signal!(
    sat::ConstantDopplerSatellite{S},
    sample_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
) where S <: AbstractGNSSSystem
    resize!(sat.carrier_code, num_samples)
    resize!(sat.signal, num_samples)
    carrier = fpcarrier!(
        sat.carrier_code,
        intermediate_frequency + get_carrier_doppler(sat),
        sample_frequency,
        get_carrier_phase_2pi(sat),
        bits = Val(7)
    )
    carrier_code = multiply_with_code!(
        carrier,
        S,
        get_code_frequency(S) + get_code_doppler(sat),
        sample_frequency,
        get_code_phase(sat),
        get_prn(sat)
    )
    sat.signal .= carrier_code .* get_amplitude(sat, n0) ./ Int16(1) << 7
end

function multiply_with_code!(
    carrier_code,
    ::Type{S},
    code_frequency,
    sample_frequency,
    start_code_phase,
    prn::Integer
) where S <: AbstractGNSSSystem
    fp = 64 - min_bits_for_code_length(S) - 1
    fp_delta = floor(Int, code_frequency * 1 << fp / sample_frequency)
    fp_start_code_phase = floor(Int, start_code_phase * 1 << fp)
    fp_total_code_length = floor(
        Int,
        get_code_length(S) * get_secondary_code_length(S) * 1 << fp
    )
    fp_code_phase = fp_start_code_phase - fp_delta
    @inbounds for i = 1:length(carrier_code)
        fp_code_phase += fp_delta
        fp_code_phase -= (fp_code_phase >= fp_total_code_length) * fp_total_code_length
        code_index = fp_code_phase >> fp
        code = get_code_unsafe(S, code_index, prn)
        carrier_code[i] *= code
    end
    carrier_code
end


"""
$(SIGNATURES)

Propagates the satellite state over time
"""
function propagate(
    sat::ConstantDopplerSatellite{S, T, D, E, C},
    num_samples,
    intermediate_frequency,
    sample_frequency,
    rng
) where {S, T, D, E, C}
    carrier_phase = (intermediate_frequency + get_carrier_doppler(sat)) * num_samples /
        sample_frequency + get_carrier_phase_2pi(sat)
    modded_carrier_phase = mod(carrier_phase + 0.5, 1) - 0.5
    code_phase = (get_code_frequency(S) + get_code_doppler(sat)) * num_samples /
        sample_frequency + get_code_phase(sat)
    modded_code_phase = mod(code_phase, get_code_length(S) * get_secondary_code_length(S))
    Δt = num_samples / sample_frequency
    exists = propagate(sat.exists, Δt, rng)
    doa = propagate(sat.doa, Δt, rng)
    ConstantDopplerSatellite{S, T, D, E, C}(
        sat.prn,
        sat.carrier_doppler,
        modded_carrier_phase,
        modded_code_phase,
        get_carrier_to_noise_density_ratio(sat),
        exists,
        doa,
        sat.carrier_code,
        sat.signal
    )
end

@inline get_amplitude(sat::AbstractSatellite, n0) = calc_amplitude_from_cn0(sat.cn0, n0)
@inline get_carrier_doppler(sat::AbstractSatellite) = sat.carrier_doppler
@inline get_code_doppler(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem =
    sat.carrier_doppler * get_code_center_frequency_ratio(S)
@inline get_carrier_phase(sat::AbstractSatellite) = 2π * sat.carrier_phase
@inline get_carrier_phase_2pi(sat::AbstractSatellite) = sat.carrier_phase
@inline get_code_phase(sat::AbstractSatellite) = sat.code_phase
@inline get_prn(sat::AbstractSatellite) = sat.prn
@inline get_gnss_system(sat::AbstractSatellite{S}) where S <: AbstractGNSSSystem = S
@inline get_carrier_to_noise_density_ratio(sat::AbstractSatellite) = sat.cn0

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
