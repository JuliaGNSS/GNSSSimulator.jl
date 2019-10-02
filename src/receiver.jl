abstract type AbstractReceiver end

struct Receiver{
        G <: Union{Number, SMatrix, AbstractGainPhaseMismCrosstalk},
        A <: Union{RotXYZ{<:Real}, AbstractAttitude}
    } <: AbstractReceiver
    sample_frequency::typeof(1.0Hz)
    intermediate_frequency::typeof(1.0Hz)
    gain_phase_mism_crosstalk::G
    attitude::A
    noise_std::Float64
end

function Receiver(
    sample_frequency;
    intermediate_frequency = 0.0Hz,
    gain_phase_mism_crosstalk = 1.0,
    attitude = RotXYZ(0,0,0),
    n0 = 1/Hz,
    noise_std = sqrt(n0 * sample_frequency)
    )
    Receiver(sample_frequency, intermediate_frequency, gain_phase_mism_crosstalk, attitude, noise_std)
end

function propagate(receiver::Receiver, Δt, rng)
    gain_phase_mism_crosstalk = propagate(receiver.gain_phase_mism_crosstalk, Δt, rng)
    attitude = propagate(receiver.attitude, Δt, rng)
    Receiver(receiver.sample_frequency, receiver.intermediate_frequency, gain_phase_mism_crosstalk, attitude, receiver.noise_std)
end

function get_gain_phase_mism_crosstalk(receiver::Receiver)
    get_gain_phase_mism_crosstalk(receiver.gain_phase_mism_crosstalk)
end

function get_attitude(receiver::Receiver)
    get_attitude(receiver.attitude)
end

@inline function get_sample_frequency(receiver::Receiver)
    receiver.sample_frequency
end

@inline function get_intermediate_frequency(receiver::Receiver)
    receiver.intermediate_frequency
end

@inline get_noise_std(receiver::Receiver) = receiver.noise_std
