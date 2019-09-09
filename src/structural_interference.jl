abstract type AbstractStructuralInterference <: AbstractEmitter end

struct ConstantDopplerStructuralInterference{
    S <: ConstantDopplerSatellite
} <: AbstractStructuralInterference
    sat::S
end

function ConstantDopplerStructuralInterference(
        sat::ConstantDopplerSatellite,
        signal_amplification;
        added_carrier_doppler = NaN*Hz,
        added_carrier_phase = NaN,
        added_code_phase = NaN,
        exists = true,
        doa = SVector(0,0,1),
        added_relative_velocity = 0.0m/s,
        added_signal_path = 0.0m
    )
    if isnan(added_carrier_doppler)
        carrier_doppler = get_carrier_doppler(sat) + added_relative_velocity / SPEED_OF_LIGHT * sat.system.center_freq
    else
        carrier_doppler = get_carrier_doppler(sat) + added_carrier_doppler
    end
    if isnan(added_carrier_phase)
        carrier_phase = mod2pi(get_carrier_phase(sat) + calc_carrier_phase(added_signal_path, sat.system.center_freq + carrier_doppler))
    else
        carrier_phase = mod2pi(get_carrier_phase(sat) + added_carrier_phase)
    end
    if isnan(added_code_phase)
        code_doppler = carrier_doppler * sat.system.code_freq / sat.system.center_freq
        code_phase = mod(sat.code_phase + calc_code_phase(added_signal_path, sat.system.code_freq + code_doppler, sat.system.code_length), sat.system.code_length)
    else
        code_phase = mod(sat.code_phase + added_code_phase, sat.system.code_length)
    end
    amplitude = get_amplitude(sat) * sqrt(uconvertp(NoUnits, signal_amplification))
    sat = ConstantDopplerSatellite(get_prn(sat), get_system(sat), carrier_doppler, carrier_phase, code_phase, amplitude, exists, doa)
    ConstantDopplerStructuralInterference(sat)
end

function propagate(si::ConstantDopplerStructuralInterference, Δt)
    ConstantDopplerStructuralInterference(propagate(si.sat, Δt))
end

function get_signal(si::AbstractStructuralInterference, attitude, get_steer_vec)
    get_signal(si.sat, attitude, get_steer_vec)
end

get_system(si::ConstantDopplerStructuralInterference) = get_system(si.sat)
get_doa(si::ConstantDopplerStructuralInterference) = get_doa(si.sat)
get_existence(si::ConstantDopplerStructuralInterference) = get_existence(si.sat)
get_carrier_doppler(si::ConstantDopplerStructuralInterference) = get_carrier_doppler(si.sat)
get_code_doppler(si::ConstantDopplerStructuralInterference) = get_code_doppler(si.sat)
get_carrier_phase(si::ConstantDopplerStructuralInterference) = get_carrier_phase(si.sat)
get_code_phase(si::ConstantDopplerStructuralInterference) = get_code_phase(si.sat)
get_prn(si::ConstantDopplerStructuralInterference) = get_prn(si.sat)
get_amplitude(si::ConstantDopplerStructuralInterference) = get_amplitude(si.sat)
