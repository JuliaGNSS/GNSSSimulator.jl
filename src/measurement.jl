"""
($SIGNATURES)

Initialize the _sim_measurement function, takes a Vector of sats and jammers `emitter`, a `gnss_system` structs an `attitude`, a steering vector function `get_steer_vec`
    and the frequencies `sample_freq` and `interm_freq`.
"""
function init_sim_measurement(emitters::Vector{T}, gnss_system::GNSSSystem, attitude::A, get_steer_vec, sample_freq, interm_freq) where T <: AbstractEmitter where A <: Union{RotXYZ, AbstractDynamicAttitude} # Die Schreibweise where A <: Union{RotXYZ, AbstractDynamicAttitude} muss möglicherweise korrigiert werden
    sim_emitter_signals = map(emitter -> init_sim_emitter_signal(emitter, gnss_system, sample_freq, interm_freq), emitters)
    init_time = 0.0
    num_samples -> _sim_measurement(num_samples, init_time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq)
end

"""
($SIGNATURES)

Return a new version of itself and a calculated `signal` with `num_samples` and the `internal_state`.
Simulates all emitters and the correlated signals of all of them.
"""
function _sim_measurement(num_samples, time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq)
    next_time = time + num_samples / sample_freq
    return_values = map(sim_emitter_signals, emitters) do sim_emitter_signal, emitter
        rotated_doa = current_attitude(time, attitude) * CartesianFromSpherical()(current_enu_doa(time, emitter.enu_doa))
        next_sim_emitter_signal, signal, internal_state = sim_emitter_signal(num_samples)
        steered_emitter_signal = get_steer_vec(rotated_doa) .* transpose(signal)
        steered_emitter_signal, next_sim_emitter_signal, internal_state
    end
    steered_emitter_signals = map(x -> x[1], return_values)
    next_sim_emitter_signals = map(x -> x[2], return_values)
    internal_states = map(x -> x[3], return_values)
    num_ants = length(get_steer_vec([0,0,1]))
    signal = sum(steered_emitter_signals) #+ gen_noise(num_ants, num_samples)
    num_samples -> _sim_measurement(num_samples, next_time, next_sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq), signal, internal_states
end
