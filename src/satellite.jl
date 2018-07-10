function init_sim_emitter_signal(sat::Satellite, gnss_system::GNSSSystem, sample_freq, interm_freq)
    sat_user_distance = init_sat_distance(sat.distance_from_earth_center, sat.enu_doa)
    init_doppler = doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, gnss_system.f₀)
    init_carrier_phase = carrier_phase(sat_user_distance, gnss_system.f₀ + init_doppler)
    init_code_phase = code_phase(sat_user_distance, gnss_system.code_f₀ + init_doppler * gnss_system.code_f₀ / gnss_system.f₀, gnss_system.code_length)
    init_amplitude = calc_amplitude_from_cn0(sat.CN0, sample_freq) # Assumes sample_freq == bandwidth
    samples -> _sim_sat_signal(samples, gnss_system, init_code_phase, init_carrier_phase, init_doppler, sample_freq, interm_freq, init_amplitude, sat.svid)
end

function _sim_sat_signal(samples, gnss_system, code_phase, carrier_phase, doppler, sample_freq, interm_freq, amplitude, svid)
    code_freq_with_doppler = gnss_system.code_f₀ + doppler * gnss_system.code_f₀ / gnss_system.f₀
    carrier_freq_with_doppler = interm_freq + doppler
    next_code_phase = gnss_system.calc_next_code_phase(samples, code_freq_with_doppler, code_phase, sample_freq)
    sampled_code = gnss_system.gen_sampled_code(1:samples, code_freq_with_doppler, code_phase, sample_freq, svid)
    next_carrier_phase = get_carrier_phase(samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    sampled_carrier = gen_carrier(1:samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    signal = sampled_code .* sampled_carrier .* amplitude
    samples -> _sim_sat_signal(samples, gnss_system, next_code_phase, next_carrier_phase, doppler, sample_freq, amplitude, svid), signal, code_phase, doppler
end
