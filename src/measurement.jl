"""
$(SIGNATURES)

Simulates post correlation measurement. It depends on:
'sat_channels' is a struct containing the channel number, 
    the complex signal after correlation, its DOA and a 
    boolean field which is '1' if the signal exists for both
    the satellite and the inference signal;
'attitudes' is a struct containing a 3×3 rotation matrix of 
    the antenna's current rotation which can be either static
    or time-dependent and optionally noisy;
'gain_phase_mism_and_crosstalk' provides the gain and phase 
    mismatch and crosstalk matrix over time 't';
'get_steer_vec' provides the steering vectors for the DOA and 
    attitude;
'noise_power' specifies the noise in [dB].
"""

function sim_post_corr_measurement(
    sat_channels,
    attitudes,
    gain_phase_mism_and_crosstalk,
    get_steer_vec,
    noise_power)

    t -> begin
        curr_attitude = sim_attitude(t, attitudes)
        𝐂 = gain_phase_mism_and_crosstalk(t) # (num_ants)x(num_ants) matrix

        sat_channel_states = map(sat_channels) do sat_channel
            sat_exists = sim_existence(t, sat_channel.exists)
            sat_doa = sim_doa(t, sat_channel.enu_doa)
            sat_signal = sim_pseudo_post_corr_signal(t, sat_channel.signal)
            interf_exists = sim_existence(t, sat_channel.interf_exists)
            interf_doa = sim_doa(t, sat_channel.interf_enu_doa)
            interf_signal = sim_pseudo_post_corr_signal(t, sat_channel.interf_signal)
            SatelliteChannelState(sat_doa, sat_signal, sat_exists, interf_doa, interf_signal, interf_exists)
        end

        measurement_wo_crosstalk = reduce(hcat, 
            map(sat_channel_states) do sat_channel_state
                sat_steering_vector = get_steer_vec(curr_attitude * sat_channel_state.doa)
                interf_steering_vector = get_steer_vec(curr_attitude * sat_channel_state.interf_doa)
                num_ants = length(sat_steering_vector)
                noise = sim_noise(noise_power, num_ants)
                sat_steering_vector .* sat_channel_state.signal .* sat_channel_state.exists .+ 
                    interf_steering_vector .* sat_channel_state.interf_signal .* sat_channel_state.interf_exists .+ 
                    noise
            end
        )

        measurement = 𝐂 * measurement_wo_crosstalk

        measurement, InternalStates(sat_channel_states, curr_attitude, 𝐂)
end
end