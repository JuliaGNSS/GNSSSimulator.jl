struct TemporalData{T}
    data::T
    sample_freq::Float64
end

function init_gen_existing_sats(
    existing_sats_over_time,
    existing_sats_sample_freq
    )
    t -> begin
        index = floor(Int, t * existing_sats_sample_freq) + 1
        existing_sats_over_time[:, index]
    end
end

function init_gen_signal_ampl_and_phase(
        max_num_sats,
        SNR_dB,
        init_signal_ampl_std,
        init_signal_phase_std,
        signal_ampl_over_time_std,
        signal_phase_over_time_std
    )

    init_signal_phase = randn(max_num_sats) * init_signal_phase_std
    init_signal_ampl = ones(max_num_sats) + randn(max_num_sats) * init_signal_ampl_std
    amplitude = 10^(SNR_dB / 20)
    (t, existing_sats) -> begin
        signal_phase = init_signal_phase + randn(max_num_sats) * signal_phase_over_time_std
        signal_ampl = init_signal_ampl + randn(max_num_sats) * signal_ampl_over_time_std
        signal_ampl[existing_sats] .* cis.(signal_phase[existing_sats]) * amplitude
    end
end

function init_gen_noise(num_ants = 1)
    (t, existing_sats) -> begin
        num_sats = sum(existing_sats)
        complex.(randn(num_ants, num_sats), randn(num_ants, num_sats)) / sqrt(2)
    end
end
