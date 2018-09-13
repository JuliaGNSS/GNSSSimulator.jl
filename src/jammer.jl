"""
$(SIGNATURES)

Simulates a jammer signal with jammer parameters `jammer`, GNSS system `gnss_system`, sample frequency `sample_freq` and intermediate frequency `interm_freq`.
Returns a function which is dependent on the number of samples `num_samples`.
"""
function init_sim_emitter_signal(jammer::CWJammer, gnss_system::S, sample_freq, interm_freq) where S <: AbstractGNSSSystem
    init_doppler = doppler(jammer.relative_velocity, gnss_system.center_freq)
    init_carrier_phase = 0.0
    init_amplitude = uconvertrp(NoUnits, jammer.JNR) * sample_freq / 1Hz
    num_samples -> _sim_cw_jammer_signal(num_samples, gnss_system, init_carrier_phase, init_doppler, sample_freq, init_amplitude, interm_freq)
end

function _sim_cw_jammer_signal(num_samples, gnss_system, carrier_phase, doppler, sample_freq, amplitude, interm_freq)
    carrier_freq_with_doppler = interm_freq + doppler
    next_carrier_phase = get_carrier_phase(num_samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    sampled_carrier = gen_carrier(1:num_samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    signal = sampled_carrier .* amplitude
    num_samples -> _sim_cw_jammer_signal(num_samples, gnss_system, next_carrier_phase, doppler, sample_freq, amplitude, interm_freq), signal, EmitterInternalStates(doppler, carrier_phase, NaN)
end
