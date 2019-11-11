abstract type AbstractJammer <: AbstractEmitter end

struct CWJammerPhase{T}
    carrier::T
end

struct CWJammerPhaseWrap
    carrier::Int
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

struct NoiseJammerPhaseWrap
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

function CWJammer(
    id,
    amplitude;
    doppler = 0.0Hz,
    phase = 0.0,
    exists = true,
    doa = SVector(0, 0, 1)
)
    phase = phase / 2π
    CWJammer(id, float(doppler), float(phase), float(amplitude), exists, doa)
end

function NoiseJammer(id, amplitude, exists = true, doa = SVector(0, 0, 1))
    NoiseJammer(id, amplitude, exists, doa)
end

function propagate(
    jammer::CWJammer,
    intermediate_frequency,
    Δt,
    rng
)
    phase_delta = (intermediate_frequency + get_carrier_doppler(jammer)) * Δt
    phase = mod(get_carrier_phase_2pi(jammer) + phase_delta + 0.5, 1) - 0.5
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    CWJammer(
        jammer.id,
        get_carrier_doppler(jammer),
        phase,
        get_amplitude(jammer),
        exists,
        doa
    )
end

function calc_phase(
    jammer::CWJammer,
    index,
    intermediate_frequency,
    sample_frequency
)
    phase = (intermediate_frequency + get_carrier_doppler(jammer)) * index /
        sample_frequency + get_carrier_phase_2pi(jammer)
    CWJammerPhase(phase)
end

function init_phase_wrap(jammer::CWJammer)
    CWJammerPhaseWrap(0)
end

function update_phase_wrap(
    phase_wrap::CWJammerPhaseWrap,
    phase::CWJammerPhase,
    jammer::CWJammer,
)
    wrap = phase_wrap.carrier + ((phase.carrier - phase_wrap.carrier) >= 0.5) * 1
    CWJammerPhaseWrap(wrap)
end

Base.@propagate_inbounds function get_signal(
    phase::CWJammerPhase,
    phase_wrap::CWJammerPhaseWrap,
    jammer::CWJammer,
    steer_vec::S,
    rng
) where {
    T <: Union{Float32, Float64},
    S <: Union{SVector{N, Complex{T}}, Complex{T}, T} where N
}
    temp = get_existence(jammer) * T(get_amplitude(jammer)) *
        GNSSSignals.cis_vfast(T(2π * phase.carrier - phase_wrap.carrier))
    steer_vec * temp
end

function propagate(jammer::NoiseJammer, intermediate_frequency, Δt, rng)
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    NoiseJammer(jammer.id, get_amplitude(jammer), exists, doa)
end

function calc_phase(
    jammer::NoiseJammer,
    index,
    intermediate_frequency,
    sample_frequency
)
    NoiseJammerPhase()
end

function init_phase_wrap(jammer::NoiseJammer)
    NoiseJammerPhaseWrap()
end

function update_phase_wrap(
    phase_wrap::NoiseJammerPhaseWrap,
    phase::NoiseJammerPhase,
    jammer::NoiseJammer,
)
    phase_wrap
end

Base.@propagate_inbounds function get_signal(
    phase::NoiseJammerPhase,
    phase_wrap::NoiseJammerPhaseWrap,
    jammer::NoiseJammer,
    steer_vec::S,
    rng
) where {
    T <: Union{Float32, Float64},
    S <: Union{SVector{N, Complex{T}}, Complex{T}, T} where N
}
    temp = get_existence(jammer) * T(get_amplitude(jammer)) * randn(rng, Complex{T})
    steer_vec * temp
end

@inline get_carrier_phase(jammer::CWJammer) = 2π * jammer.phase
@inline get_carrier_phase_2pi(jammer::CWJammer) = jammer.phase
@inline get_carrier_doppler(jammer::CWJammer) = jammer.doppler
@inline get_id(jammer::AbstractJammer) = jammer.id
