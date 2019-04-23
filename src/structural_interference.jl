"""
$(SIGNATURES)

Simulates a structural interference for one satellite signal with satellite parameters `sat`, GNSS system `gnss_system`, sample frequency `sample_freq` and intermediate frequency `interm_freq`.
Returns a function which is dependent on the number of samples `num_samples`.
"""
function init_sim_emitter_signal(interference::StructuralInterference, gnss_system::S, sample_freq, interm_freq) where S <: AbstractGNSSSystem
    init_sat_user_distance = calc_init_sat_user_distance(interference.sat.distance_from_earth_center, interference.sat.enu_doa)
    init_doppler = calc_init_doppler(interference.sat.distance_from_earth_center, interference.sat.enu_doa, interference.sat.velocity, gnss_system.center_freq) + interference.added_relative_velocity / SPEED_OF_LIGHT * gnss_system.center_freq
    init_carrier_phase = calc_carrier_phase(init_sat_user_distance + interference.added_signal_path, gnss_system.center_freq + init_doppler)
    init_code_phase = calc_code_phase(init_sat_user_distance + interference.added_signal_path, gnss_system.code_freq + init_doppler * gnss_system.code_freq / gnss_system.center_freq, gnss_system.code_length)
    init_amplitude = calc_amplitude_from_cn0(interference.sat.CN0 + interference.signal_amplification, 1/1Hz) # Assumes sample_freq == bandwidth
    num_samples -> _sim_sat_signal(num_samples, gnss_system, init_code_phase, init_carrier_phase, init_doppler, sample_freq, interm_freq, init_amplitude, interference.sat.prn), EmitterInternalStates(init_doppler, init_carrier_phase, init_code_phase)
end
