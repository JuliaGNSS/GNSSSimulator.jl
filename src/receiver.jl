abstract type AbstractReceiver end

struct Receiver{
        G <: Union{Float64, SMatrix, AbstractGainPhaseMismCrosstalk},
        A <: Union{RotXYZ{Float64}, AbstractAttitude},
        Q <: Noise
    } <: AbstractReceiver
    gain_phase_mism_crosstalk::G
    attitude::A
    noise::Q
end

function Receiver(gain_phase_mism_crosstalk::Float64, attitude, noise_std::Float64)
    noise = Noise(randn(ComplexF64) * noise_std, noise_std)
    Receiver(gain_phase_mism_crosstalk, attitude, noise)
end

function Receiver(gain_phase_mism_crosstalk, attitude, noise_std::Float64)
    gpmc = get_gain_phase_mism_crosstalk(gain_phase_mism_crosstalk)
    noise = Noise(randn(SVector{size(gpmc, 1), ComplexF64}) * noise_std, noise_std)
    Receiver(gain_phase_mism_crosstalk, attitude, noise)
end

function propagate(receiver::Receiver, Δt)
    gain_phase_mism_crosstalk = propagate(receiver.gain_phase_mism_crosstalk, Δt)
    attitude = propagate(receiver.attitude, Δt)
    noise = propagate(receiver.noise, Δt)
    Receiver(gain_phase_mism_crosstalk, attitude, noise)
end

function get_gain_phase_mism_crosstalk(receiver::Receiver)
    get_gain_phase_mism_crosstalk(receiver.gain_phase_mism_crosstalk)
end

function get_attitude(receiver::Receiver)
    get_attitude(receiver.attitude)
end

function get_noise(receiver::Receiver)
    get_noise(receiver.noise)
end
