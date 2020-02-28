abstract type AbstractReceiver{T} end

struct Receiver{
        T <: AbstractFloat,
        G <: Union{T, SMatrix{N, N, Complex{T}}, AbstractGainPhaseMismCrosstalk{T}} where N,
        A <: Union{RotXYZ{<:Real}, AbstractAttitude},
        C <: Unitful.Quantity{T, Unitful.𝐓}
    } <: AbstractReceiver{T}
    sample_frequency::typeof(1.0Hz)
    intermediate_frequency::typeof(1.0Hz)
    gain_phase_mism_crosstalk::G
    attitude::A
    n0::C
end

function Receiver(
    sample_frequency;
    intermediate_frequency = 0.0Hz,
    gain_phase_mism_crosstalk::G = 1.0,
    attitude = RotXYZ(0,0,0),
    n0::Unitful.Quantity{<:Real, Unitful.𝐓} = 1/Hz,
) where {
    N,
    T <: AbstractFloat,
    G <: Union{T, SMatrix{N, N, Complex{T}}, AbstractGainPhaseMismCrosstalk{T}}
}
    Receiver(
        sample_frequency,
        intermediate_frequency,
        gain_phase_mism_crosstalk,
        attitude,
        T(n0.val)/Hz
    )
end

function propagate(receiver::Receiver, num_samples, rng)
    Δt = num_samples / get_sample_frequency(receiver)
    gain_phase_mism_crosstalk = propagate(receiver.gain_phase_mism_crosstalk, Δt, rng)
    attitude = propagate(receiver.attitude, Δt, rng)
    Receiver(
        receiver.sample_frequency,
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

@inline function get_sample_frequency(receiver::Receiver)
    receiver.sample_frequency
end

@inline function get_intermediate_frequency(receiver::Receiver)
    receiver.intermediate_frequency
end

@inline get_noise_density(receiver::AbstractReceiver) = receiver.n0
@inline function get_noise_std(receiver::AbstractReceiver{T}) where T <: AbstractFloat
    T(sqrt(get_noise_density(receiver) * get_sample_frequency(receiver)))
end
