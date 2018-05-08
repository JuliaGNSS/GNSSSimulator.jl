function init_measurement(
        attitudes,
        existing_sats,
        sat_doa_carts,
        get_steer_vec;
        SNR_dB = 15,
        init_phase_mism_betw_ant_std = π / 4,
        phase_mism_over_time_std = 1,
        init_gain_mism_betw_ant_std = 0.1,
        gain_mism_over_time_std = 0.01,
        init_crosstalk_to_direct_power_dB = -15,
        init_crosstalk_ampl_std = 0.05,
        init_crosstalk_phase_std = 2 * π,
        crosstalk_ampl_over_time_std = 0.001,
        crosstalk_phase_over_time_std = 0.1,
        attitude_over_time_std = 1,
        init_signal_ampl_std = 0.05,
        init_signal_phase_std = 2 * π,
        signal_ampl_over_time_std = 0.1,
        signal_phase_over_time_std = 0.5
    )

    num_ants = size(get_steer_vec(Spherical(SVector(0.0,0.0,1.0))), 1)
    max_num_sats = size(existing_sats, 1)

    gen_gain_and_phase_mism_and_crosstalk = init_gen_gain_and_phase_mism_and_crosstalk(
        num_ants,
        init_phase_mism_betw_ant_std,
        phase_mism_over_time_std,
        init_gain_mism_betw_ant_std,
        gain_mism_over_time_std,
        init_crosstalk_to_direct_power_dB,
        init_crosstalk_ampl_std,
        init_crosstalk_phase_std,
        crosstalk_ampl_over_time_std,
        crosstalk_phase_over_time_std)

    gen_signal_ampl_and_phase = init_gen_signal_ampl_and_phase(
        max_num_sats,
        SNR_dB,
        init_signal_ampl_std,
        init_signal_phase_std,
        signal_ampl_over_time_std,
        signal_phase_over_time_std)

    gen_attitude = init_gen_attitude(
        attitudes.data,
        attitudes.sample_freq,
        attitude_over_time_std)

    gen_steering_vectors = init_gen_steering_vectors(get_steer_vec)

    gen_doas = init_gen_doas(
        sat_doa_carts.data,
        sat_doa_carts.sample_freq)

    gen_existing_sats = init_gen_existing_sats(
        existing_sats.data,
        existing_sats.sample_freq)

    gen_noise = init_gen_noise(num_ants)

    t -> begin
        existing_sats = gen_existing_sats(t)
        attitude = gen_attitude(t)
        doas = gen_doas(t, existing_sats)
        𝐀 = gen_steering_vectors(t, attitude, doas)
        𝐂 = gen_gain_and_phase_mism_and_crosstalk(t)
        𝐬 = gen_signal_ampl_and_phase(t, existing_sats)
        𝐍 = gen_noise(t, existing_sats)
        𝐘 = 𝐂 * (𝐀 .* 𝐬.' + 𝐍)
        𝐘, attitude, doas, 𝐀, 𝐬, 𝐂, existing_sats
    end
end

function init_gen_gain_and_phase_mism_and_crosstalk(
        num_ants,
        init_phase_mism_betw_ant_std,
        phase_mism_over_time_std,
        init_gain_mism_betw_ant_std,
        gain_mism_over_time_std,
        init_crosstalk_to_direct_power_dB,
        init_crosstalk_ampl_std,
        init_crosstalk_phase_std,
        crosstalk_ampl_over_time_std,
        crosstalk_phase_over_time_std
    )

    init_phase_mism = randn(num_ants) * init_phase_mism_betw_ant_std * π / 180
    init_gain_mism = randn(num_ants) * init_gain_mism_betw_ant_std

    init_crosstalk_ampl = (ones(num_ants, num_ants) + randn(num_ants, num_ants) * init_crosstalk_ampl_std) .*
        (ones(num_ants, num_ants) - eye(num_ants)) * 10^(init_crosstalk_to_direct_power_dB / 10)
    init_crosstalk_phase = randn(num_ants, num_ants) * init_crosstalk_phase_std

    t -> begin
        phase_mism = init_phase_mism + randn(num_ants) * phase_mism_over_time_std * π / 180
        gain_mism = init_gain_mism + randn(num_ants) * gain_mism_over_time_std
        gain_and_phase_mism = gain_mism .* cis.(phase_mism)
        crosstalk_phase = init_crosstalk_phase + randn(num_ants, num_ants) * crosstalk_phase_over_time_std * π / 180
        crosstalk_ampl = init_crosstalk_ampl + randn(num_ants, num_ants) * crosstalk_ampl_over_time_std .* (ones(num_ants, num_ants) - eye(num_ants))
        crosstalk = crosstalk_ampl .* cis.(crosstalk_phase)

        normalize_gain_and_phase_mism_and_crosstalk(diagm(gain_and_phase_mism) + crosstalk)
    end
end

function normalize_gain_and_phase_mism_and_crosstalk(gain_and_phase_mism_and_crosstalk)
    gain_and_phase_mism_and_crosstalk_norm = mapslices(norm, gain_and_phase_mism_and_crosstalk, 1)
    gain_and_phase_mism_and_crosstalk ./ gain_and_phase_mism_and_crosstalk_norm
end

function init_gen_doas(
        sats_doa_over_time,
        sats_doa_sample_freq
    )
    (t, existing_sats) -> begin
        index = floor(Int, t * sats_doa_sample_freq) + 1
        sats_doa_over_time[:,existing_sats,index]
    end
end

function init_gen_steering_vectors(
        get_steer_vec
    )
    (t, attitude, doas) -> begin
        mapslices(doa -> get_steer_vec(Spherical(attitude * doa)), doas, 1)
    end
end

function init_gen_attitude(
        attitude_over_time,
        attitude_sample_freq,
        attitude_over_time_std
    )
    t -> begin
        index = floor(Int, t * attitude_sample_freq) + 1
        RotXYZ(attitude_over_time[:,index]...) * RotXYZ(randn(SVector{3}) * attitude_over_time_std * π / 180...)
    end
end
