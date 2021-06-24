abstract type AbstractJammer{T} <: AbstractEmitter{T} end

struct CWJammer{
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Union{Unitful.Gain, AbstractAmplitude},
} <: AbstractJammer{T}
    id::Int
    doppler::typeof(1.0Hz)
    phase::Float64
    jnr::C
    exists::E
    doa::D
    signal::StructArray{Complex{T},1,NamedTuple{(:re, :im),Tuple{Array{T,1},Array{T,1}}},Int}
end

struct NoiseJammer{
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Union{Unitful.Gain, AbstractAmplitude},
} <: AbstractJammer{T}
    id::Int
    jnr::C
    exists::E
    doa::D
    signal::StructArray{Complex{T},1,NamedTuple{(:re, :im),Tuple{Array{T,1},Array{T,1}}},Int}
end

function CWJammer(
    ::Type{T},
    id::Integer,
    jnr::J;
    doppler = 0.0Hz,
    phase = 0.0,
    exists::E = true,
    doa::D = SVector(0, 0, 1)
) where {
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    J <: Union{Unitful.Gain{Unitful.Decibel, :?, <:Real}, AbstractAmplitude}
}
    phase = phase / 2π
    signal = StructArray{Complex{T}}(undef, 0)
    CWJammer{T, D, E, J}(
        id,
        float(doppler),
        float(phase),
        jnr,
        exists,
        doa,
        signal
    )
end

function CWJammer(
    id::Integer,
    jnr::Union{Unitful.Gain{Unitful.Decibel, :?, <:Real}, AbstractAmplitude};
    doppler = 0.0Hz,
    phase = 0.0,
    exists = true,
    doa = SVector(0, 0, 1)
)
    CWJammer(
        Float32,
        id,
        jnr,
        doppler = doppler,
        phase = phase,
        exists = exists,
        doa = doa
    )
end

function gen_signal!(
    jammer::CWJammer{T},
    sampling_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
) where T
    resize!(jammer.signal, num_samples)
    gen_carrier!(
        jammer.signal,
        intermediate_frequency + get_carrier_doppler(jammer),
        sampling_frequency,
        get_carrier_phase_2pi(jammer),
        T(get_amplitude(jammer.jnr, n0, sampling_frequency))
    )
end

function NoiseJammer(
    ::Type{T},
    id,
    jnr::J;
    exists::E = true,
    doa::D = SVector(0, 0, 1)
) where {
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    J <: Union{Unitful.Gain{Unitful.Decibel, :?, <:Real}, AbstractAmplitude}
}
    signal = StructArray{Complex{T}}(undef, 0)
    NoiseJammer{T, D, E, J}(
        id,
        jnr,
        exists,
        doa,
        signal
    )
end

function NoiseJammer(
    id,
    jnr;
    exists = true,
    doa = SVector(0, 0, 1)
)
    NoiseJammer(
        Float32,
        id,
        jnr,
        exists = exists,
        doa = doa
    )
end

function gen_signal!(
    jammer::NoiseJammer{T},
    sample_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
) where T
    resize!(jammer.signal, num_samples)
    signal = randn!(rng, jammer.signal)
    jammer.signal .*= T(get_amplitude(jammer.jnr, n0, sample_frequency))
    jammer.signal
end

function propagate(
    jammer::CWJammer,
    num_samples,
    intermediate_frequency,
    sample_frequency,
    n0,
    rng
)
    carrier_phase = (intermediate_frequency + get_carrier_doppler(jammer)) * num_samples /
        sample_frequency + get_carrier_phase_2pi(jammer)
    modded_carrier_phase = mod(carrier_phase + 0.5, 1) - 0.5
    Δt = num_samples / sample_frequency
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    jnr = propagate(jammer.jnr, n0, sample_frequency, Δt, rng)
    CWJammer(
        jammer.id,
        get_carrier_doppler(jammer),
        modded_carrier_phase,
        jnr,
        exists,
        doa,
        jammer.signal
    )
end

function propagate(
    jammer::NoiseJammer,
    num_samples,
    intermediate_frequency,
    sample_frequency,
    n0,
    rng
)
    Δt = num_samples / sample_frequency
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    jnr = propagate(jammer.jnr, n0, sample_frequency, Δt, rng)
    NoiseJammer(jammer.id, jnr, exists, doa, jammer.signal)
end

@inline get_carrier_phase(jammer::CWJammer) = 2π * jammer.phase
@inline get_carrier_phase_2pi(jammer::CWJammer) = jammer.phase
@inline get_carrier_doppler(jammer::CWJammer) = jammer.doppler
@inline get_id(jammer::AbstractJammer) = jammer.id
@inline get_jammer_to_noise_ratio(jammer::AbstractJammer) = get_jammer_to_noise_ratio(jammer.jnr)
