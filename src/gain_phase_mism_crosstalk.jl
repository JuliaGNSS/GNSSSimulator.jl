abstract type AbstractGainPhaseMismCrosstalk end

"""
(a * t^e + b) / (t^e + 1)
"""
struct AsymptoticGainPhaseMismCrosstalk{
        N,
        E <: Union{Float64, Int64},
        T <: SMatrix{N, N, ComplexF64}
    } <: AbstractGainPhaseMismCrosstalk
    a::SVector{N, Float64}
    b::SVector{N, Float64}
    e::E
    C::T
    t::typeof(1.0s)
end

function AsymptoticGainPhaseMismCrosstalk(
    a::SVector{N, Float64},
    b::SVector{N, Float64};
    e = 2,
    C = SMatrix{N,N,ComplexF64}(I),
    t = 0.0s
) where N
    AsymptoticGainPhaseMismCrosstalk(a, b, e, C, t)
end

function propagate(gpmc::AsymptoticGainPhaseMismCrosstalk, Δt, rng)
    t = gpmc.t + Δt
    AsymptoticGainPhaseMismCrosstalk(gpmc.a, gpmc.b, gpmc.e, gpmc.C, t)
end

function get_gain_phase_mism_crosstalk(gpmc::AsymptoticGainPhaseMismCrosstalk)
    phase_misms = cis.(
        (gpmc.a * (gpmc.t / s) ^ gpmc.e + gpmc.b) ./ ((gpmc.t / s) ^ gpmc.e + 1)
    )
    phase_misms .* gpmc.C
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
