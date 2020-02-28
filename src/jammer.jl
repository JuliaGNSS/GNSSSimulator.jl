abstract type AbstractJammer{T} <: AbstractEmitter{T} end

struct CWJammer{
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Unitful.Gain{Unitful.Decibel, :?, T}
} <: AbstractJammer{T}
    id::Int
    doppler::typeof(1.0Hz)
    phase::Float64
    jnr::C
    exists::E
    doa::D
    carrier::StructArray{Complex{Int16},1,NamedTuple{(:re, :im),Tuple{Array{Int16,1},Array{Int16,1}}},Int}
    signal::StructArray{Complex{T},1,NamedTuple{(:re, :im),Tuple{Array{T,1},Array{T,1}}},Int}
end

struct NoiseJammer{
    T <: AbstractFloat,
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
    C <: Unitful.Gain{Unitful.Decibel, :?, T}
} <: AbstractJammer{T}
    id::Int
    jnr::C
    exists::E
    doa::D
    signal::StructArray{Complex{T},1,NamedTuple{(:re, :im),Tuple{Array{T,1},Array{T,1}}},Int}
end

function CWJammer(
    id::Integer,
    jnr::Unitful.Gain{Unitful.Decibel, :?, <:Real};
    doppler = 0.0Hz,
    phase = 0.0,
    exists::E = true,
    doa::D = SVector(0, 0, 1)
) where {
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence},
}
    phase = phase / 2π
    carrier = StructArray{Complex{Int16}}(undef, 0)
    signal = StructArray{Complex{eltype(float(jnr.val))}}(undef, 0)
    CWJammer{eltype(float(jnr.val)), D, E, eltype(float(jnr))}(
        id,
        float(doppler),
        float(phase),
        float(jnr),
        exists,
        doa,
        carrier,
        signal
    )
end

function gen_signal!(
    jammer::CWJammer,
    sample_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
)
    resize!(jammer.carrier, num_samples)
    resize!(jammer.signal, num_samples)
    carrier = fpcarrier!(
        jammer.carrier,
        intermediate_frequency + get_carrier_doppler(jammer),
        sample_frequency,
        get_carrier_phase_2pi(jammer),
        bits = Val(7)
    )
    jammer.signal .= carrier .* get_amplitude(jammer, n0, sample_frequency) ./ 1 << 7
    jammer.signal
end

function NoiseJammer(
    id,
    jnr::Unitful.Gain{Unitful.Decibel, :?, <:Real};
    exists::E = true,
    doa::D = SVector(0, 0, 1)
) where {
    D <: Union{SVector{3}, AbstractDOA},
    E <: Union{Bool, AbstractExistence}
}
    signal = StructArray{Complex{eltype(float(jnr.val))}}(undef, 0)
    NoiseJammer{eltype(float(jnr.val)), D, E, eltype(float(jnr))}(
        id,
        float(jnr),
        exists,
        doa,
        signal
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
    jammer.signal .= signal .* get_amplitude(jammer, n0, sample_frequency)
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
        jammer.carrier,
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

@inline function get_amplitude(jammer::AbstractJammer, n0, sample_frequency)
    sqrt(uconvertp(NoUnits, jammer.jnr) * n0 * sample_frequency)
end
@inline get_carrier_phase(jammer::CWJammer) = 2π * jammer.phase
@inline get_carrier_phase_2pi(jammer::CWJammer) = jammer.phase
@inline get_carrier_doppler(jammer::CWJammer) = jammer.doppler
@inline get_id(jammer::AbstractJammer) = jammer.id
@inline get_jammer_to_noise_ratio(jammer::AbstractJammer) = jammer.jnr
