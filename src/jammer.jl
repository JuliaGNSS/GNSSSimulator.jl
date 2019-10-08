abstract type AbstractJammer <: AbstractEmitter end

struct CWJammerPhase{T}
    carrier::T
end

struct CWJammer{
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
} <: AbstractJammer
    id::Int
    doppler::typeof(1.0Hz)
    phase::Float64
    amplitude::Float64
    exists::E
    doa::D
end

struct NoiseJammerPhase
end

struct NoiseJammer{
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
} <: AbstractJammer
    id::Int
    amplitude::Float64
    exists::E
    doa::D
end

function CWJammer(id, amplitude; doppler = 0.0Hz, phase = 0.0, exists = true, doa = SVector(0, 0, 1))
    CWJammer(id, doppler, phase, amplitude, exists, doa)
end

function NoiseJammer(id, amplitude, exists = true, doa = SVector(0, 0, 1))
    NoiseJammer(id, amplitude, exists, doa)
end

function propagate(jammer::CWJammer, intermediate_frequency, Δt, rng) where T <: Union{Float32, Float64}
    phase_delta = 2π * upreferred((intermediate_frequency + get_carrier_doppler(jammer)) * Δt)
    phase = mod2pi(get_carrier_phase(jammer) + phase_delta + π) - π
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    CWJammer(jammer.id, get_carrier_doppler(jammer), phase, get_amplitude(jammer), exists, doa)
end

function fast_propagate(phase::CWJammerPhase, jammer::CWJammer, intermediate_frequency, Δt)
    carrier_delta = 2π * upreferred((intermediate_frequency + get_carrier_doppler(jammer)) * Δt)
    carrier_phase = phase.carrier + carrier_delta
    carrier_phase = carrier_phase - (carrier_phase > π) * 2π
    CWJammerPhase(carrier_phase)
end

Base.@propagate_inbounds function get_signal(phase::CWJammerPhase, jammer::CWJammer, steer_vec::S, rng) where {T <: Union{Float32, Float64}, S <: Union{SVector{N, Complex{T}}, Complex{T}, T} where N}
    temp = get_existence(jammer) * T(get_amplitude(jammer)) * GNSSSignals.cis_vfast(T(phase.carrier))
    steer_vec * temp
end

function propagate(jammer::NoiseJammer, intermediate_frequency, Δt, rng)
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    NoiseJammer(jammer.id, get_amplitude(jammer), exists, doa)
end

function fast_propagate(phase::NoiseJammerPhase, jammer::NoiseJammer, intermediate_frequency, Δt)
    phase
end

Base.@propagate_inbounds function get_signal(phase::NoiseJammerPhase, jammer::NoiseJammer, steer_vec::S, rng) where {T <: Union{Float32, Float64}, S <: Union{SVector{N, Complex{T}}, Complex{T}, T} where N}
    temp = get_existence(jammer) * T(get_amplitude(jammer)) * randn(rng, Complex{T})
    steer_vec * temp
end

@inline get_carrier_phase(jammer::CWJammer) = jammer.phase
@inline get_carrier_doppler(jammer::CWJammer) = jammer.doppler
@inline get_phase(jammer::CWJammer) = CWJammerPhase(get_carrier_phase(jammer))
@inline get_phase(jammer::NoiseJammer) = NoiseJammerPhase()
@inline get_id(jammer::AbstractJammer) = jammer.id
