function init_sim_measurement(emitters::Vector{T}, gnss_system::GNSSSystem, attitude::A, get_steer_vec, sample_freq, interm_freq) where T <: AbstractEmitter where A <: Union{RotXYZ, AbstractDynamicAttitude} # Die Schreibweise where A <: Union{RotXYZ, AbstractDynamicAttitude} muss möglicherweise korrigiert werden
    sim_emitter_signals = map(emitter -> init_sim_emitter_signal(emitter, gnss_system, sample_freq, interm_freq), emitters)
    init_time = 0.0
    samples -> _sim_measurement(samples, init_time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq)
end

function _sim_measurement(samples, time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq)
    next_time = time + samples / sample_freq
    steered_emitter_signals = map(sim_emitter_signals, emitters) do sim_emitter_signal, emitter
        rotated_doa = current_attitude(time, attitude) * current_enu_doa(time, emitter.enu_doa)
        get_steer_vec(rotated_doa) .* sim_emitter_signal(samples)
    end
    num_ants = length(get_steer_vec([0,0,1]))
    signal = sum(steered_emitter_signals) + gen_noise(num_ants, samples)
    samples -> _sim_measurement(samples, next_time, sim_emitter_signals, emitters, attitude, get_steer_vec, sample_freq), signal
end
