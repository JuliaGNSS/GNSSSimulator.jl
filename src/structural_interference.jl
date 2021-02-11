abstract type AbstractStructuralInterference{T} <: AbstractEmitter{T} end

struct ConstantDopplerStructuralInterference{
    T <: AbstractFloat,
    CS <: ConstantDopplerSatellite{<: AbstractGNSS, T}
} <: AbstractStructuralInterference{T}
    sat::CS
end

function ConstantDopplerStructuralInterference(
    sat::ConstantDopplerSatellite{S, T},
    signal_amplification::Unitful.Gain{Unitful.Decibel, :?, <:Real};
    added_carrier_doppler = NaN*Hz,
    added_carrier_phase = NaN,
    added_code_phase = NaN,
    exists::E = true,
    doa::D = SVector(0,0,1),
    added_relative_velocity = 0.0m/s,
    added_signal_path = 0.0m
) where {
    S <: AbstractGNSS,
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
}
    system = sat.system
    if isnan(added_carrier_doppler)
        carrier_doppler = get_carrier_doppler(sat) + added_relative_velocity /
            SPEED_OF_LIGHT * get_center_frequency(system)
    else
        carrier_doppler = get_carrier_doppler(sat) + added_carrier_doppler
    end
    if isnan(added_carrier_phase)
        carrier_phase = mod(get_carrier_phase_2pi(sat) +
            calc_carrier_phase(
                added_signal_path,
                get_center_frequency(system) + carrier_doppler
            ) + 0.5,
            1
        ) - 0.5
    else
        added_carrier_phase = added_carrier_phase / 2π
        carrier_phase = mod(get_carrier_phase_2pi(sat) + added_carrier_phase + 0.5, 1) - 0.5
    end
    if isnan(added_code_phase)
        code_doppler = carrier_doppler * get_code_center_frequency_ratio(system)
        code_phase = mod(
            get_code_phase(sat) +
            calc_code_phase(
                added_signal_path,
                get_code_frequency(system) + code_doppler, get_code_length(system)
            ),
            get_code_length(system)
        )
    else
        code_phase = mod(get_code_phase(sat) + added_code_phase, get_code_length(system))
    end
    cn0 = get_carrier_to_noise_density_ratio(sat) + signal_amplification
    code = Vector{Int8}(undef, 0)
    signal = StructArray{Complex{T}}(undef, 0)
    sat = ConstantDopplerSatellite{S, T, D, E, eltype(float(cn0))}(
        system,
        get_prn(sat),
        carrier_doppler,
        carrier_phase,
        code_phase,
        cn0,
        exists,
        doa,
        code,
        signal
    )
    ConstantDopplerStructuralInterference(sat)
end

function gen_signal!(
    si::ConstantDopplerStructuralInterference,
    sample_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
)
    gen_signal!(si.sat, sample_frequency, intermediate_frequency, n0, num_samples, rng)
end

function propagate(
    si::ConstantDopplerStructuralInterference,
    num_samples,
    intermediate_frequency,
    sample_frequency,
    rng
)
    ConstantDopplerStructuralInterference(
        propagate(si.sat, num_samples, intermediate_frequency, sample_frequency, rng)
    )
end

@inline get_doa(si::ConstantDopplerStructuralInterference) = get_doa(si.sat)
@inline get_existence(si::ConstantDopplerStructuralInterference) = get_existence(si.sat)
@inline get_carrier_doppler(si::ConstantDopplerStructuralInterference) =
    get_carrier_doppler(si.sat)
@inline get_code_doppler(si::ConstantDopplerStructuralInterference) =
    get_code_doppler(si.sat)
@inline get_carrier_phase(si::ConstantDopplerStructuralInterference) =
    get_carrier_phase(si.sat)
@inline get_code_phase(si::ConstantDopplerStructuralInterference) = get_code_phase(si.sat)
@inline get_prn(si::ConstantDopplerStructuralInterference) = get_prn(si.sat)
@inline get_amplitude(si::ConstantDopplerStructuralInterference, n0) =
    get_amplitude(si.sat, n0)
@inline get_gnss_system(si::ConstantDopplerStructuralInterference) = get_gnss_system(si.sat)
@inline get_carrier_to_noise_density_ratio(si::ConstantDopplerStructuralInterference) =
    get_carrier_to_noise_density_ratio(si.sat)
