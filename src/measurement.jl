"""
$(SIGNATURES)

Simulates the measurement based on a Array of `emitters` with a subtype of `AbstractEmitter`. `emitters` may also be a scalar.
The GNSS System is specified by `gnss_system`, the attitude of the antenna by `attitude`, the sample frequency by `sample_freq` and 
the intermediate frequency by `interm_freq`. The steering vectors for each DOA are specified by the function `get_steer_vec` which
depends on the DOA. The added white noise can be disabled by `add_noise`.
"""
function init_sim_measurement(emitters::Vector{T}, gnss_system::GNSSSystem, attitude::A, get_steer_vec, sample_freq, interm_freq, add_noise::Bool = true) where T <: AbstractEmitter where A <: Union{RotXYZ, AbstractDynamicAttitude} # Die Schreibweise where A <: Union{RotXYZ, AbstractDynamicAttitude} muss möglicherweise korrigiert werden
    sim_emitter_signals = map(emitter -> init_sim_emitter_signal(emitter, gnss_system, sample_freq, interm_freq), emitters)
    init_time = 0.0
    num_samples -> _sim_measurement(num_samples, init_time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq, add_noise)
end

function _sim_measurement(num_samples, time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq, add_noise::Bool)
    next_time = time + num_samples / sample_freq
    return_values = map(sim_emitter_signals, emitters) do sim_emitter_signal, emitter
        cart_enu_doa = CartesianFromSpherical()(current_enu_doa(time, emitter.enu_doa))
        rotated_doa = current_attitude(time, attitude) * cart_enu_doa
        next_sim_emitter_signal, signal, internal_state = sim_emitter_signal(num_samples)
        steered_emitter_signal = get_steer_vec(rotated_doa) .* transpose(signal)
        steered_emitter_signal, next_sim_emitter_signal, internal_state
    end
    steered_emitter_signals = map(x -> x[1], return_values)
    next_sim_emitter_signals = map(x -> x[2], return_values)
    internal_states = map(x -> x[3], return_values)
    num_ants = length(get_steer_vec([0,0,1]))
    signal = sum(steered_emitter_signals) .+ (add_noise ? gen_noise(num_ants, num_samples) : 0)
    num_samples -> _sim_measurement(num_samples, next_time, next_sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq, add_noise), signal, internal_states
end
