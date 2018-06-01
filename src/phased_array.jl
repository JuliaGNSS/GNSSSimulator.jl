struct InternalStates
        doas::Matrix{Float64}
        existing_sats::Vector{Float64}
        attitude::RotXYZ
        gain_phase_mism_crosstalk::Matrix{Complex{Float64}}
        steering_vectors::Matrix{Complex{Float64}}
        signal::Vector{Complex{Float64}}
end

function init_measurement(
        attitudes,
        existing_sats,
        sat_doa_carts,
        get_steer_vec;
        SNR_dB = 15,
        init_phase_mism_betw_ant_std = π / 12,
        phase_mism_over_time_std = π / 180,
        init_gain_mism_betw_ant_std = 0.01,
        gain_mism_over_time_std = 0.001,
        init_crosstalk_to_direct_power_dB = -15,
        init_crosstalk_ampl_std = 0.05,
        init_crosstalk_phase_std = 2 * π,
        crosstalk_ampl_over_time_std = 0.001,
        crosstalk_phase_over_time_std = 0.1 * π / 180,
        yaw_over_time_std = 0.1 * π / 180,
        pitch_over_time_std = 0.01 * π / 180,
        roll_over_time_std = 0.01 * π / 180,
        init_signal_ampl_std = 0.05,
        init_signal_phase_std = 2 * π,
        signal_ampl_over_time_std = 0.01,
        signal_phase_over_time_std = 0.5 * π / 180
    )

    num_ants = size(get_steer_vec(Spherical(SVector(0.0,0.0,1.0))), 1)
    max_num_sats = size(existing_sats.data, 1)

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
        init_signal_ampl_std,
        init_signal_phase_std,
        signal_ampl_over_time_std,
        signal_phase_over_time_std)

    gen_attitude = init_gen_attitude(
        attitudes.data,
        attitudes.sample_freq,
        yaw_over_time_std,
        pitch_over_time_std,
        roll_over_time_std)

    gen_steering_vectors = init_gen_steering_vectors(get_steer_vec)

    gen_doas = init_gen_doas(
        sat_doa_carts.data,
        sat_doa_carts.sample_freq)

    gen_existing_sats = init_gen_existing_sats(
        existing_sats.data,
        existing_sats.sample_freq)

    gen_noise = init_gen_noise(-SNR_dB, num_ants)

    t -> _measurement(t, gen_existing_sats, gen_attitude, gen_doas, gen_steering_vectors, gen_gain_and_phase_mism_and_crosstalk, gen_signal_ampl_and_phase, gen_noise)
end

function _measurement(t, gen_existing_sats, gen_attitude, gen_doas, gen_steering_vectors, gen_gain_and_phase_mism_and_crosstalk, gen_signal_ampl_and_phase, gen_noise)
    existing_sats = gen_existing_sats(t)
    attitude = gen_attitude(t)
    doas = gen_doas(t, existing_sats)
    𝐀 = gen_steering_vectors(t, attitude, doas)
    𝐂 = gen_gain_and_phase_mism_and_crosstalk(t)
    𝐬 = gen_signal_ampl_and_phase(t, existing_sats)
    𝐍 = gen_noise(t, existing_sats)
    𝐘::Array{Complex{Float64}, 2} = 𝐂 * (𝐀 .* 𝐬.' + 𝐍)
    internal_states = InternalStates(doas, existing_sats, attitude, 𝐂, 𝐀, 𝐬)
    𝐘, internal_states
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

    init_phase_mism = randn(num_ants) * init_phase_mism_betw_ant_std
    init_gain_mism = ones(num_ants) + randn(num_ants) * init_gain_mism_betw_ant_std

    init_crosstalk_ampl = (ones(num_ants, num_ants) + randn(num_ants, num_ants) * init_crosstalk_ampl_std) .*
        (ones(num_ants, num_ants) - eye(num_ants)) * 10^(init_crosstalk_to_direct_power_dB / 10)
    init_crosstalk_phase = randn(num_ants, num_ants) * init_crosstalk_phase_std

    t -> begin
        phase_mism = init_phase_mism + randn(num_ants) * phase_mism_over_time_std
        gain_mism = init_gain_mism + randn(num_ants) * gain_mism_over_time_std
        gain_and_phase_mism = gain_mism .* cis.(phase_mism)
        crosstalk_phase = init_crosstalk_phase + randn(num_ants, num_ants) * crosstalk_phase_over_time_std
        crosstalk_ampl = init_crosstalk_ampl + randn(num_ants, num_ants) * crosstalk_ampl_over_time_std .* (ones(num_ants, num_ants) - eye(num_ants))
        crosstalk = crosstalk_ampl .* cis.(crosstalk_phase)

        normalize_gain_and_phase_mism_and_crosstalk(diagm(gain_and_phase_mism) + crosstalk)
    end
end

function normalize_gain_and_phase_mism_and_crosstalk(gain_and_phase_mism_and_crosstalk)
    gain_and_phase_mism_and_crosstalk_norm = map(norm, julienne(gain_and_phase_mism_and_crosstalk, (:,*)))'
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
        hcat(map(get_steer_vec, julienne(attitude * doas, (:,*)))...)
    end
end

function init_gen_attitude(
        attitude_over_time,
        attitude_sample_freq,
        yaw_over_time_std,
        pitch_over_time_std,
        roll_over_time_std
    )
    t -> begin
        index = floor(Int, t * attitude_sample_freq) + 1
        yaw_noise = randn() * yaw_over_time_std
        pitch_noise = randn() * pitch_over_time_std
        roll_noise = randn() * roll_over_time_std
        curr_att = attitude_over_time[:,index]
        RotXYZ(curr_att[1] + roll_noise, curr_att[2] + pitch_noise, curr_att[3] + yaw_noise)
    end
end
