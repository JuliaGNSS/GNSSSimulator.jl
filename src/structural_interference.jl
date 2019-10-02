abstract type AbstractStructuralInterference <: AbstractEmitter end

struct ConstantDopplerStructuralInterference{
    S <: ConstantDopplerSatellite
} <: AbstractStructuralInterference
    sat::S
end

function ConstantDopplerStructuralInterference(
        sat::ConstantDopplerSatellite{S, D, E},
        signal_amplification;
        added_carrier_doppler = NaN*Hz,
        added_carrier_phase = NaN,
        added_code_phase = NaN,
        amplitude = NaN,
        exists = true,
        doa = SVector(0,0,1),
        added_relative_velocity = 0.0m/s,
        added_signal_path = 0.0m
    )  where {S <: AbstractGNSSSystem,  D <: Union{SVector{3}, AbstractDOA}, E <: Union{Bool, AbstractExistence}}
    if isnan(added_carrier_doppler)
        carrier_doppler = get_carrier_doppler(sat) + added_relative_velocity / SPEED_OF_LIGHT * get_center_frequency(S)
    else
        carrier_doppler = get_carrier_doppler(sat) + added_carrier_doppler
    end
    if isnan(added_carrier_phase)
        carrier_phase = mod2pi(get_carrier_phase(sat) + calc_carrier_phase(added_signal_path, get_center_frequency(S) + carrier_doppler))
    else
        carrier_phase = mod2pi(get_carrier_phase(sat) + added_carrier_phase)
    end
    if isnan(added_code_phase)
        code_doppler = carrier_doppler * get_code_center_frequency_ratio(S)
        code_phase = mod(get_code_phase(sat) + calc_code_phase(added_signal_path, get_code_frequency(S) + code_doppler, get_code_length(S)), get_code_length(S))
    else
        code_phase = mod(get_code_phase(sat) + added_code_phase, get_code_length(S))
    end
    if isnan(amplitude)
        amplitude = get_amplitude(sat) * sqrt(uconvertp(NoUnits, signal_amplification))
    end
    sat = ConstantDopplerSatellite{S, D, E}(get_prn(sat), carrier_doppler, carrier_phase, code_phase, amplitude, exists, doa)
    ConstantDopplerStructuralInterference(sat)
end

function fast_propagate(phase::SatellitePhase, si::ConstantDopplerStructuralInterference, intermediate_frequency, Δt)
    fast_propagate(phase, si.sat, intermediate_frequency, Δt)
end

Base.@propagate_inbounds function get_signal(phase::SatellitePhase, si::ConstantDopplerStructuralInterference, steer_vec, rng)
    get_signal(phase, si.sat, steer_vec, rng)
end

function propagate(si::ConstantDopplerStructuralInterference, phase, Δt, rng)
    ConstantDopplerStructuralInterference(propagate(si.sat, phase, Δt, rng))
end

@inline get_doa(si::ConstantDopplerStructuralInterference) = get_doa(si.sat)
@inline get_existence(si::ConstantDopplerStructuralInterference) = get_existence(si.sat)
@inline get_carrier_doppler(si::ConstantDopplerStructuralInterference) = get_carrier_doppler(si.sat)
@inline get_code_doppler(si::ConstantDopplerStructuralInterference) = get_code_doppler(si.sat)
@inline get_carrier_phase(si::ConstantDopplerStructuralInterference) = get_carrier_phase(si.sat)
@inline get_code_phase(si::ConstantDopplerStructuralInterference) = get_code_phase(si.sat)
@inline get_prn(si::ConstantDopplerStructuralInterference) = get_prn(si.sat)
@inline get_amplitude(si::ConstantDopplerStructuralInterference) = get_amplitude(si.sat)
@inline get_phase(si::ConstantDopplerStructuralInterference) = SatellitePhase(get_carrier_phase(si), get_code_phase(si))
