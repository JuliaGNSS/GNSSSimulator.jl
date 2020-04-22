abstract type AbstractGainPhaseMismCrosstalk{T} end

"""
(a * t^e + b) / (t^e + 1)
"""
struct AsymptoticGainPhaseMismCrosstalk{
        N,
        T <: AbstractFloat,
        E <: Union{Float64, Int64},
        CT <: SMatrix{N, N, Complex{T}}
    } <: AbstractGainPhaseMismCrosstalk{T}
    a::SVector{N, T}
    b::SVector{N, T}
    e::E
    C::CT
    t::typeof(1.0s)
end

function AsymptoticGainPhaseMismCrosstalk(
    a::SVector{N, T},
    b::SVector{N, T};
    e = 2,
    C = SMatrix{N,N,Complex{T}}(I),
    t = 0.0s
) where {N, T <: AbstractFloat}
    AsymptoticGainPhaseMismCrosstalk(a, b, e, C, t)
end

function propagate(gpmc::AsymptoticGainPhaseMismCrosstalk, Δt, rng)
    t = gpmc.t + Δt
    AsymptoticGainPhaseMismCrosstalk(gpmc.a, gpmc.b, gpmc.e, gpmc.C, t)
end

function get_gain_phase_mism_crosstalk(gpmc::AsymptoticGainPhaseMismCrosstalk{N,T}) where {N,T}
    phase_misms = cis.(
        (gpmc.a * (gpmc.t / s) ^ gpmc.e + gpmc.b) ./ ((gpmc.t / s) ^ gpmc.e + 1)
    )
    SMatrix{N,N,Complex{T}}(phase_misms .* gpmc.C)
end

@inline function propagate(gain_phase_mism_crosstalk, Δt, rng)
    gain_phase_mism_crosstalk
end

@inline function get_gain_phase_mism_crosstalk(gain_phase_mism_crosstalk)
    gain_phase_mism_crosstalk
end

"""
$(SIGNATURES)

Normalizes the gain and phase mismatch and crosstalk function so that the norm over columns is 1.

"""
function normalize_gain_phase_mism_and_crosstalk(gain_and_phase_mism_and_crosstalk)
    gain_and_phase_mism_and_crosstalk_norm = map(
        norm,
        eachcol(gain_and_phase_mism_and_crosstalk)
    )'
    gain_and_phase_mism_and_crosstalk ./ gain_and_phase_mism_and_crosstalk_norm
end
