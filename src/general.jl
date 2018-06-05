"""
$(SIGNATURES)

Simulates the direction of arrivals over time `t`. This is the default configuration from Lothar Kurz.

"""
function sim_doas()
    doas = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
           -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
            0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]
    (t, existing_sats) -> doas[:,existing_sats]
end

"""
$(SIGNATURES)

Simulates the direction over time based on the data `doas_over_time` with the sample frequency `sample_freq`.
The first dimension of the data should contain the cartesian unit vectors. The second dimension should be the time.

"""
function sim_doas(doas_over_time, sample_freq)
    (t, existing_sats) -> begin
        index = floor(Int, t * sample_freq) + 1
        doas_over_time[:,existing_sats,index]
    end
end

"""
$(SIGNATURES)

Simulates a static satellite existence over time `t` based on the boolean array `existing_sats`. 

"""
function sim_existing_sats(existing_sats)
    t -> existing_sats
end

"""
$(SIGNATURES)

Simulates a varying satellite existence over time `t` based on the data `existing_sats_over_time`
and sample frequency `sample_freq`. The first dimension should hold existence of satellites and
the second dimension the time.

"""
function sim_existing_sats(existing_sats_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        existing_sats_over_time[:, index]
    end
end

"""
$(SIGNATURES)

Simulates a pseudo post correlation signal over time `t` for given `existing_sats` at that time 
instance with the power of `signal_power_dB`. Pseudo means that no GNSS data is included.

"""
function sim_pseudo_post_corr_signal(num_sats, signal_power_dB, init_phase_var_between_signals = π)
    amplitude = 10^(signal_power_dB / 20)
    init_signal_phase = randn(num_sats) * sqrt(init_phase_var_between_signals)
    signal = amplitude .* cis.(init_signal_phase)
    (t, existing_sats) -> begin
        signal[existing_sats]
    end
end

"""
$(SIGNATURES)

Simulates a varying pseudo post correlation signal over time `t` for given `existing_sats` at that time 
instance with the power of `signal_power_dB`. Pseudo means that no GNSS data is included.

"""
function sim_pseudo_post_corr_signal(num_sats, signal_power_dB, init_phase_var_between_signals, ampl_var, phase_var)
    pseudo_post_corr_signal = sim_pseudo_post_corr_signal(num_sats, signal_power_dB, init_phase_var_between_signals)
    ampl_std = sqrt(ampl_var)
    phase_std = sqrt(phase_var)
    (t, existing_sats) -> begin
        curr_num_sats = sum(existing_sats)
        signal = pseudo_post_corr_signal(t, existing_sats)
        (abs.(signal) .+ randn(curr_num_sats) .* ampl_std) .* cis.(angle.(signal) .+ randn(curr_num_sats) .* phase_std)
    end
end

"""
$(SIGNATURES)

Simulates a static attitude over time `t`.

"""
function sim_attitude(init_yaw, init_pitch, init_roll)
    attitude = RotXYZ(init_roll, init_pitch, init_yaw)
    t -> attitude
end

"""
$(SIGNATURES)

Simulates a varying attitude over time `t` given by the data `attitude_over_time` with sample frequency
`sample_freq`. The first dimension should hold the attitude angles in the order of roll, pitch and yaw.
The second dimension should hold the time.

"""
function sim_attitude(attitude_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        RotXYZ(attitude_over_time[1, index], attitude_over_time[2, index], attitude_over_time[3, index])
    end
end

"""
$(SIGNATURES)

Simulates a varying attitude over time `t` given by the data `attitude_over_time` with sample frequency
`sample_freq` and variances given by `yaw_over_time_var`, `pitch_over_time_var` and `roll_over_time_var`. 
The first dimension should hold the attitude angles in the order of roll, pitch and yaw. The second 
dimension should hold the time.

"""
function sim_attitude(attitude_over_time, sample_freq, yaw_over_time_var, pitch_over_time_var, roll_over_time_var)
    yaw_over_time_std = sqrt(yaw_over_time_var)
    pitch_over_time_std = sqrt(pitch_over_time_var)
    roll_over_time_std = sqrt(roll_over_time_var)
    t -> begin
        yaw_noise = randn() * yaw_over_time_std
        pitch_noise = randn() * pitch_over_time_std
        roll_noise = randn() * roll_over_time_std
        index = floor(Int, t * sample_freq) + 1
        curr_att = attitude_over_time[:,index]
        RotXYZ(curr_att[1] + roll_noise, curr_att[2] + pitch_noise, curr_att[3] + yaw_noise)
    end
end

"""
$(SIGNATURES)

Simulates noise over time `t` for the `existing_sats` at that time instance with the power of
`noise_power_dB` and number of antennas `num_ants`.

"""
function sim_noise(noise_power_dB, num_ants = 1)
    amplitude = 10^(noise_power_dB / 20)
    (t, existing_sats) -> begin
        num_sats = sum(existing_sats)
        complex.(randn(num_ants, num_sats), randn(num_ants, num_sats)) / sqrt(2) * amplitude
    end
end
