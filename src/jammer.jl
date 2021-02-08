abstract type AbstractJammer{T} <: AbstractEmitter{T} end

struct CWJammer{
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Unitful.Gain
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
    C <: Unitful.Gain
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
    jnr::Unitful.Gain{Unitful.Decibel, :?, <:Real};
    doppler = 0.0Hz,
    phase = 0.0,
    exists::E = true,
    doa::D = SVector(0, 0, 1)
) where {
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
}
    phase = phase / 2π
    signal = StructArray{Complex{T}}(undef, 0)
    CWJammer{T, D, E, eltype(float(jnr))}(
        id,
        float(doppler),
        float(phase),
        float(jnr),
        exists,
        doa,
        signal
    )
end

function CWJammer(
    id::Integer,
    jnr::Unitful.Gain{Unitful.Decibel, :?, <:Real};
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
    jammer::CWJammer,
    sampling_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
)
    resize!(jammer.signal, num_samples)
    gen_carrier!(
        jammer.signal,
        intermediate_frequency + get_carrier_doppler(jammer),
        sampling_frequency,
        get_carrier_phase_2pi(jammer),
        get_amplitude(jammer, n0, sampling_frequency)
    )
end

function NoiseJammer(
    ::Type{T},
    id,
    jnr::Unitful.Gain{Unitful.Decibel, :?, <:Real};
    exists::E = true,
    doa::D = SVector(0, 0, 1)
) where {
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
}
    signal = StructArray{Complex{T}}(undef, 0)
    NoiseJammer{T, D, E, eltype(float(jnr))}(
        id,
        float(jnr),
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
    jammer::NoiseJammer,
    sample_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
)
    resize!(jammer.signal, num_samples)
    signal = randn!(rng, jammer.signal)
    jammer.signal .*= get_amplitude(jammer, n0, sample_frequency)
    jammer.signal
end

function propagate(
    jammer::CWJammer,
    num_samples,
    intermediate_frequency,
    sample_frequency,
    rng
)
    carrier_phase = (intermediate_frequency + get_carrier_doppler(jammer)) * num_samples /
        sample_frequency + get_carrier_phase_2pi(jammer)
    modded_carrier_phase = mod(carrier_phase + 0.5, 1) - 0.5
    Δt = num_samples / sample_frequency
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    CWJammer(
        jammer.id,
        get_carrier_doppler(jammer),
        modded_carrier_phase,
        get_jammer_to_noise_ratio(jammer),
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
    rng
)
    Δt = num_samples / sample_frequency
    exists = propagate(jammer.exists, Δt, rng)
    doa = propagate(jammer.doa, Δt, rng)
    NoiseJammer(jammer.id, get_jammer_to_noise_ratio(jammer), exists, doa, jammer.signal)
end

@inline function get_amplitude(jammer::AbstractJammer{T}, n0, sample_frequency) where T
    T(sqrt(uconvertp(NoUnits, jammer.jnr) * n0 * sample_frequency))
end
@inline get_carrier_phase(jammer::CWJammer) = 2π * jammer.phase
@inline get_carrier_phase_2pi(jammer::CWJammer) = jammer.phase
@inline get_carrier_doppler(jammer::CWJammer) = jammer.doppler
@inline get_id(jammer::AbstractJammer) = jammer.id
@inline get_jammer_to_noise_ratio(jammer::AbstractJammer) = jammer.jnr
