"""
$(SIGNATURES)

Calculate amplitude of a signal based on the carrier-to-noise-density-ratio (CN0) `cn0` in [dB-Hz] and the frequency bandwidth `bandwidth` in [Hz], assumes noise power to be 1.
"""
function calc_amplitude_from_cn0(cn0, n0)
    sqrt(linear(cn0) * n0)
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

Simulates a satellite signal with satellite parameters `sat`, GNSS system `gnss_system`, sample frequency `sample_freq` and intermediate frequency `interm_freq`.
Returns a function which is dependent on the number of samples `num_samples`.
"""
function init_sim_emitter_signal(sat::Satellite, gnss_system::S, sample_freq, interm_freq) where S <: AbstractGNSSSystem
    init_sat_user_distance = calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    init_doppler = calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, gnss_system.center_freq)
    init_carrier_phase = calc_carrier_phase(init_sat_user_distance, gnss_system.center_freq + init_doppler)
    init_code_phase = calc_code_phase(init_sat_user_distance, gnss_system.code_freq + init_doppler * gnss_system.code_freq / gnss_system.center_freq, gnss_system.code_length)
    init_amplitude = calc_amplitude_from_cn0(sat.CN0, 1/1Hz) # Assumes sample_freq == bandwidth
    num_samples -> _sim_sat_signal(num_samples, gnss_system, init_code_phase, init_carrier_phase, init_doppler, sample_freq, interm_freq, init_amplitude, sat.prn)
end

function _sim_sat_signal(num_samples, gnss_system, code_phase, carrier_phase, doppler, sample_freq, interm_freq, amplitude, prn)
    code_freq_with_doppler = gnss_system.code_freq + doppler * gnss_system.code_freq / gnss_system.center_freq
    carrier_freq_with_doppler = interm_freq + doppler
    next_code_phase = GNSSSignals.calc_code_phase(num_samples, code_freq_with_doppler, code_phase, sample_freq, gnss_system.code_length)
    sampled_code = gen_code(gnss_system, 1:num_samples, code_freq_with_doppler, code_phase, sample_freq, prn)
    next_carrier_phase = get_carrier_phase(num_samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    sampled_carrier = gen_carrier(1:num_samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    signal = sampled_code .* sampled_carrier .* amplitude
    num_samples -> _sim_sat_signal(num_samples, gnss_system, next_code_phase, next_carrier_phase, doppler, sample_freq, interm_freq, amplitude, prn), signal, EmitterInternalStates(doppler, carrier_phase, code_phase)
end
