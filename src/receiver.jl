abstract type AbstractReceiver end

struct Receiver{
        G <: Union{Float64, SMatrix, AbstractGainPhaseMismCrosstalk},
        A <: Union{RotXYZ{Float64}, AbstractAttitude}
    } <: AbstractReceiver
    gain_phase_mism_crosstalk::G
    attitude::A
    noise_std::Float64
end

function propagate(receiver::Receiver, Δt)
    gain_phase_mism_crosstalk = propagate(receiver.gain_phase_mism_crosstalk, Δt)
    attitude = propagate(receiver.attitude, Δt)
    Receiver(gain_phase_mism_crosstalk, attitude, receiver.noise_std)
end

function get_gain_phase_mism_crosstalk(receiver::Receiver)
    get_gain_phase_mism_crosstalk(receiver.gain_phase_mism_crosstalk)
end

function get_attitude(receiver::Receiver)
    get_attitude(receiver.attitude)
end
