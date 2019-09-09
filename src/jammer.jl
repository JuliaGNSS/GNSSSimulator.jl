abstract type AbstractJammer <: AbstractEmitter end

struct CWJammer{
    T <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
} <: AbstractJammer
    id::Int
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    amplitude::Float64
    exists::E
    doa::T
end

struct NoiseJammer{
    T <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
} <: AbstractJammer
    id::Int
    noise::ComplexF64
    amplitude::Float64
    exists::E
    doa::T
end

function NoiseJammer(id, amplitude, exists = true, doa = SVector(0.0, 0.0, 1.0))
    noise = randn(ComplexF64)
    NoiseJammer(id, noise, amplitude, exists, doa)
end

function propagate(jammer::CWJammer, Δt)
    carrier_phase = 2π * jammer.carrier_doppler * Δt + jammer.carrier_phase
    carrier_phase -= (carrier_phase > 2π) * 2π
    exists = propagate(jammer.exists, Δt)
    doa = propagate(jammer.doa, Δt)
    CWJammer(jammer.id, jammer.carrier_doppler, carrier_phase, jammer.amplitude, exists, doa)
end

function propagate(jammer::NoiseJammer, Δt)
    noise = randn(ComplexF64)
    exists = propagate(jammer.exists, Δt)
    doa = propagate(jammer.doa, Δt)
    NoiseJammer(jammer.id, noise, jammer.amplitude, exists, doa)
end

get_carrier_phase(jammer::CWJammer) = jammer.carrier_phase
get_carrier_doppler(jammer::CWJammer) = jammer.carrier_doppler
get_noise(jammer::NoiseJammer) = jammer.noise

function get_signal(jammer::CWJammer, attitude, get_steer_vec)
    get_existence(jammer) * (get_steer_vec(get_doa(jammer), attitude) .* (
        cis(get_carrier_phase(jammer)) *
        get_amplitude(jammer)
    ))
end

function get_signal(jammer::NoiseJammer, attitude, get_steer_vec)
    get_existence(jammer) * (get_steer_vec(get_doa(jammer), attitude) .* (
        get_noise(jammer) *
        get_amplitude(jammer)
    ))
end
