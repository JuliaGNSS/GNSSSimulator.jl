function init_sim_emitter_signal(jammer::CWJammer, gnss_system::GNSSSystem, sample_freq, interm_freq)
    init_doppler = doppler(jammer.relative_velocity, gnss_system.f₀)
    init_carrier_phase = 0.0
    init_amplitude = calc_amplitude_from_jnr(jammer.JNR)
    samples -> _sim_cw_jammer_signal(samples, gnss_system, init_carrier_phase, init_doppler, sample_freq, init_amplitude, interm_freq)
end

function _sim_cw_jammer_signal(samples, gnss_system, carrier_phase, doppler, sample_freq, amplitude, interm_freq)
    carrier_freq_with_doppler = interm_freq + doppler
    next_carrier_phase = get_carrier_phase(samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    sampled_carrier = gen_carrier(1:samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    signal = sampled_carrier .* amplitude
    samples -> _sim_cw_jammer_signal(samples, gnss_system, next_carrier_phase, doppler, sample_freq, amplitude, interm_freq), signal
end
