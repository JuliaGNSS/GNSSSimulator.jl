"""
$(SIGNATURES)

Simulates post correlation measurement. It depends on:
`sat_channels` is a struct containing the channel number,
    the complex signal after correlation, its DOA and a
    boolean field which is '1' if the signal exists for both
    the satellite and the inference signal;
`attitudes` is a struct containing a 3×3 rotation matrix of
    the antenna's current rotation which can be either static
    or time-dependent and optionally noisy;
`gain_phase_mism_and_crosstalk` provides the gain and phase
    mismatch and crosstalk matrix over time `t`;
`get_steer_vec` provides the steering vectors for the DOA and
    attitude;
`noise_power` specifies the noise in [dB].
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

"""
$(SIGNATURES)

Simulates the measurement based on an Array of `emitters` with a subtype of `AbstractEmitter`. `emitters` may also be a scalar.
The GNSS System is specified by `gnss_system`, the attitude of the antenna by `attitude`, the sample frequency by `sample_freq` and
the intermediate frequency by `interm_freq`. The steering vectors for each DOA are specified by the function `get_steer_vec` which
depends on the DOA. The added white noise can be disabled by `add_noise`.
"""
function init_sim_measurement(emitters::Vector{T}, gnss_system::S, attitude::A, get_steer_vec, sample_freq, interm_freq, add_noise::Bool = true) where {T <: AbstractEmitter, A <: Union{RotXYZ, NoisyStaticAttitude, AbstractDynamicAttitude}, S <: AbstractGNSSSystem}
    sim_emitter_signals = map(emitter -> init_sim_emitter_signal(emitter, gnss_system, sample_freq, interm_freq), emitters)
    #sim_emitter_signals = map(emitter -> func_wrapper(init_sim_emitter_signal(emitter, gnss_system, sample_freq, interm_freq)), emitters)
    init_time = 0.0s
    num_samples -> _sim_measurement(num_samples, init_time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq, add_noise)
end

function _sim_measurement(num_samples, time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq, add_noise::Bool)
    next_time = time + num_samples / sample_freq
    return_values = map(sim_emitter_signals, emitters) do sim_emitter_signal, emitter
        cart_enu_doa = sim_doa(time, emitter.enu_doa)
        rotated_doa = sim_attitude(time, attitude) * cart_enu_doa
        exists = sim_existence(time, emitter.exists)
        next_sim_emitter_signal, emitter_signal, internal_state = sim_emitter_signal(num_samples)
        steered_emitter_signal = transpose(get_steer_vec(rotated_doa)) .* emitter_signal .* exists
        steered_emitter_signal, next_sim_emitter_signal, internal_state
    end
    steered_emitter_signals = map(x -> x[1], return_values)
    next_sim_emitter_signals = map(x -> x[2], return_values)
    internal_states = map(x -> x[3], return_values)
    num_ants = length(get_steer_vec([0,0,1]))
    signal = sum(steered_emitter_signals) .+ (add_noise ? gen_noise(num_ants, num_samples) * sqrt(sample_freq / 1Hz) : 0)
    num_samples -> _sim_measurement(num_samples, next_time, next_sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq, add_noise), signal, internal_states
end

#const func_wrapper = FunctionWrappers.FunctionWrapper{Tuple{Function, Array{Complex{Float64}, 1}, EmitterInternalStates}, Tuple{Int64}}
