struct ReceivedSignal{
    ES <: Tuple,
    R <: AbstractReceiver
}
    emitters::ES
    receiver::R
end

function get_measurement(state::ReceivedSignal, get_steer_vec = (doa, attitude) -> 1)
    attitude = get_attitude(state.receiver)
    signal = mapreduce(emitter -> get_signal(emitter, attitude, get_steer_vec), +, state.emitters)
    noise = get_noise(state.receiver)
    C = get_gain_phase_mism_crosstalk(state.receiver)
    C * (signal + noise)
end

function propagate(state::ReceivedSignal, Δt)
    emitters = map(emitter -> propagate(emitter, Δt), state.emitters)
    receiver = propagate(state.receiver, Δt)
    ReceivedSignal(emitters, receiver)
end
