struct CWJammer{T <: Union{Spherical{R}, AbstractDynamicDOA{R}} where R<:Real} <: AbstractJammer
    id::Int
    enu_doa::T
    relative_velocity::Float64
    JNR::Float64
end

"""
$(SIGNATURES)

Calculate amplitude of a signal based on the Jammer-to-Noise-Ratio (JNR) `jnr_dB` in decibel, assumes noise power to be 1.
"""
function calc_amplitude_from_jnr(jnr_dB)
    10^(jnr_dB / 20)
end

"""
$(SIGNATURES)

Simulates a jammer signal with jammer parameters `jammer`, GNSS system `gnss_system`, sample frequency `sample_freq` and intermediate frequency `interm_freq`.
Returns a function which is dependent on the number of samples `num_samples`.
"""
function init_sim_emitter_signal(jammer::CWJammer, gnss_system::GNSSSystem, sample_freq, interm_freq)
    init_doppler = doppler(jammer.relative_velocity, gnss_system.f₀)
    init_carrier_phase = 0.0
    init_amplitude = calc_amplitude_from_jnr(jammer.JNR)
    num_samples -> _sim_cw_jammer_signal(num_samples, gnss_system, init_carrier_phase, init_doppler, sample_freq, init_amplitude, interm_freq)
end

function _sim_cw_jammer_signal(num_samples, gnss_system, carrier_phase, doppler, sample_freq, amplitude, interm_freq)
    carrier_freq_with_doppler = interm_freq + doppler
    next_carrier_phase = get_carrier_phase(num_samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    sampled_carrier = gen_carrier(1:num_samples, carrier_freq_with_doppler, carrier_phase, sample_freq)
    signal = sampled_carrier .* amplitude
    num_samples -> _sim_cw_jammer_signal(num_samples, gnss_system, next_carrier_phase, doppler, sample_freq, amplitude, interm_freq), signal, EmitterInternalStates(doppler, carrier_phase, NaN)
end
