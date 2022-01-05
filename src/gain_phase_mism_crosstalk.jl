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

struct RandomWalkGainPhaseMismCrosstalk{N, G <: SMatrix, C <: SMatrix}
    gain_phase::G
    crosstalk::C
    gain_phase_vel_std::Float64
    crosstalk_vel_std::Float64
    crosstalk_ampl::Float64
end

function RandomWalkGainPhaseMismCrosstalk(N; gain_phase_vel_std, crosstalk_vel_std, crosstalk_ampl)
    gain_phase = zeros(SMatrix{2, N - 1, Float64})
    crosstalk = zeros(SMatrix{2, N^2 - N, Float64})
    RandomWalkGainPhaseMismCrosstalk{N, typeof(gain_phase), typeof(crosstalk)}(
        gain_phase,
        crosstalk,
        gain_phase_vel_std,
        crosstalk_vel_std,
        crosstalk_ampl
    )
end

@inline function propagate(gpmc::RandomWalkGainPhaseMismCrosstalk{N}, Δt, rng) where N
    P = get_process(Order(2), Δt)
    C = cholesky(get_process_covariance(Order(2), Δt)).L
    next_gain_phase = P * gpmc.gain_phase + C * randn(rng, SMatrix{2, N - 1, Float64}) * gpmc.gain_phase_vel_std
    next_crosstalk = P * gpmc.crosstalk + C * randn(rng, SMatrix{2, N^2 - N, Float64}) * gpmc.crosstalk_vel_std
    RandomWalkGainPhaseMismCrosstalk{N, typeof(next_gain_phase), typeof(next_crosstalk)}(next_gain_phase, next_crosstalk, gpmc.gain_phase_vel_std, gpmc.crosstalk_vel_std, gpmc.crosstalk_ampl)
end

@inline function get_gain_phase_mism_crosstalk(gpmc::RandomWalkGainPhaseMismCrosstalk{N}) where N
    C = MMatrix{N,N,ComplexF64}(undef)
    crosstalk_idx = 1
    for i = 1:N, j = 1:N
        if i == j == N
            C[i,i] = complex(1.0, 0.0)
        elseif i == j
            C[i,i] = cis(gpmc.gain_phase[1,i])
        else
            C[i,j] = cis(gpmc.crosstalk[1,crosstalk_idx]) * gpmc.crosstalk_ampl
            crosstalk_idx += 1
        end
    end
    SMatrix(C)
end