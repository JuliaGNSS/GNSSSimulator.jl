abstract type AbstractReceiver end

struct Receiver{
        G <: Union{TG, SMatrix{N, N, Complex{TG}}, AbstractGainPhaseMismCrosstalk{TG}} where {N, TG <: AbstractFloat},
        A <: Union{RotXYZ{<:Real}, AbstractAttitude},
        C <: Unitful.Quantity{<:Real, Unitful.𝐓}
    } <: AbstractReceiver
    sampling_frequency::typeof(1.0Hz)
    intermediate_frequency::typeof(1.0Hz)
    gain_phase_mism_crosstalk::G
    attitude::A
    n0::C
end

function Receiver(
    sampling_frequency;
    intermediate_frequency = 0.0Hz,
    gain_phase_mism_crosstalk::G = 1.0,
    attitude::A = RotXYZ(0, 0, 0),
    n0::C = 1/Hz,
) where {
    N,
    G <: Union{TG, SMatrix{N, N, Complex{TG}}, AbstractGainPhaseMismCrosstalk{TG}} where TG <: AbstractFloat,
    A <: Union{RotXYZ{<:Real}, AbstractAttitude},
    C <: Unitful.Quantity{<:Real, Unitful.𝐓}
}
    Receiver{G, A, C}(
        sampling_frequency,
        intermediate_frequency,
        gain_phase_mism_crosstalk,
        attitude,
        n0
    )
end

function propagate(receiver::Receiver, num_samples, rng)
    Δt = num_samples / get_sampling_frequency(receiver)
    gain_phase_mism_crosstalk = propagate(receiver.gain_phase_mism_crosstalk, Δt, rng)
    attitude = propagate(receiver.attitude, Δt, rng)
    Receiver(
        receiver.sampling_frequency,
        receiver.intermediate_frequency,
        gain_phase_mism_crosstalk,
        attitude,
        receiver.n0
    )
end

function get_gain_phase_mism_crosstalk(receiver::Receiver)
    get_gain_phase_mism_crosstalk(receiver.gain_phase_mism_crosstalk)
end

function get_attitude(receiver::Receiver)
    get_attitude(receiver.attitude)
end

@inline function get_sampling_frequency(receiver::Receiver)
    receiver.sampling_frequency
end

@inline function get_intermediate_frequency(receiver::Receiver)
    receiver.intermediate_frequency
end

@inline get_noise_density(receiver::AbstractReceiver) = receiver.n0
@inline function get_noise_std(receiver::AbstractReceiver)
    sqrt(get_noise_density(receiver) * get_sampling_frequency(receiver))
end
