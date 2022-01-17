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

struct RandomWalkGainPhaseMismCrosstalk{
    N,
    T <: AbstractFloat,
    P <: SMatrix,
    G <: SVector,
    CP <: SMatrix,
    CA <: SVector
} <: AbstractGainPhaseMismCrosstalk{T}
    phase_mism::P
    gain_mism::G
    crosstalk_phase::CP
    crosstalk_ampl::CA
    phase_mism_vel_std::Float64
    crosstalk_phase_vel_std::Float64
end

function RandomWalkGainPhaseMismCrosstalk(
    N,
    T = Float64;
    gain_mism_std = 0.1,
    relative_crosstalk_ampl_std = 0.1,
    phase_mism_vel_std,
    crosstalk_phase_vel_std,
    crosstalk_ampl, 
)
    gain_mism = ones(SVector{N, T}) + randn(SVector{N, T}) * gain_mism_std
    phase_mism = zeros(SMatrix{2, N, T})
    crosstalk_phase = zeros(SMatrix{2, N^2 - N, T})
    crosstalk_ampl = (ones(SVector{N^2 - N, T}) + randn(SVector{N^2 - N, T}) * relative_crosstalk_ampl_std) * crosstalk_ampl
    RandomWalkGainPhaseMismCrosstalk{N, T, typeof(phase_mism), typeof(gain_mism), typeof(crosstalk_phase), typeof(crosstalk_ampl)}(
        phase_mism,
        gain_mism,
        crosstalk_phase,
        crosstalk_ampl,
        phase_mism_vel_std,
        crosstalk_phase_vel_std,
    )
end

@inline function propagate(gpmc::RandomWalkGainPhaseMismCrosstalk{N, T, P, G, CP, CA}, Δt, rng) where {N, T, P, G, CP, CA}
    process = get_process(Order(2), upreferred(Δt / s))
    C = cholesky(get_process_covariance(Order(2), upreferred(Δt / s))).L
    next_phase_mism = process * gpmc.phase_mism + C * randn(rng, SMatrix{2, N, T}) * gpmc.phase_mism_vel_std
    next_crosstalk_phase = process * gpmc.crosstalk_phase + C * randn(rng, SMatrix{2, N^2 - N, T}) * gpmc.crosstalk_phase_vel_std
    RandomWalkGainPhaseMismCrosstalk{N, T, P, G, CP, CA}(next_phase_mism, gpmc.gain_mism, next_crosstalk_phase, gpmc.crosstalk_ampl, gpmc.phase_mism_vel_std, gpmc.crosstalk_phase_vel_std)
end

@inline function get_gain_phase_mism_crosstalk(gpmc::RandomWalkGainPhaseMismCrosstalk{N, T}) where {N, T}
    C = MMatrix{N,N,Complex{T}}(undef)
    crosstalk_idx = 1
    for i = 1:N, j = 1:N
        if i == j
            C[i,i] = cis(gpmc.phase_mism[1, i]) * gpmc.gain_mism[i]
        else
            C[i,j] = cis(gpmc.crosstalk_phase[1, crosstalk_idx]) * gpmc.crosstalk_ampl[crosstalk_idx]
            crosstalk_idx += 1
        end
    end
    SMatrix(C)
end