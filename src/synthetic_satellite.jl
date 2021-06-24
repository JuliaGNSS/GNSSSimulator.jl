abstract type AbstractSyntheticSatellite{T} <: AbstractEmitter{T} end

struct SyntheticSatellite{
    T <: AbstractFloat,
    CS <: ConstantDopplerSatellite{<: AbstractGNSS, T}
} <: AbstractSyntheticSatellite{T}
    sat::CS
end

function gen_signal!(
    synthetic_sat::SyntheticSatellite,
    sample_frequency,
    intermediate_frequency,
    n0,
    num_samples::Integer,
    rng
)
    gen_signal!(
        synthetic_sat.sat,
        sample_frequency,
        intermediate_frequency,
        n0,
        num_samples,
        rng
    )
end

function propagate(
    synthetic_sat::SyntheticSatellite,
    num_samples,
    intermediate_frequency,
    sample_frequency,
    n0,
    rng
)
    SyntheticSatellite(
        propagate(
            synthetic_sat.sat,
            num_samples,
            intermediate_frequency,
            sample_frequency,
            n0,
            rng
        )
    )
end

function get_steer_vec(
    manifold::AbstractManifold{N},
    emitter::SyntheticSatellite{T},
    attitude
) where {N, T <: AbstractFloat}
    ones(SVector{N, T})
end

function get_steer_vec(
    manifold::IdealManifold{1},
    emitter::SyntheticSatellite{T},
    attitude
) where T <: AbstractFloat
    one(T)
end


@inline get_doa(ss::SyntheticSatellite) = get_doa(ss.sat)
@inline get_existence(ss::SyntheticSatellite) = get_existence(ss.sat)
@inline get_carrier_doppler(ss::SyntheticSatellite) = get_carrier_doppler(ss.sat)
@inline get_code_doppler(ss::SyntheticSatellite) = get_code_doppler(ss.sat)
@inline get_carrier_phase(ss::SyntheticSatellite) = get_carrier_phase(ss.sat)
@inline get_code_phase(ss::SyntheticSatellite) = get_code_phase(ss.sat)
@inline get_prn(ss::SyntheticSatellite) = get_prn(ss.sat)
@inline get_amplitude(ss::SyntheticSatellite, n0) = get_amplitude(ss.sat, n0)
@inline get_gnss_system(ss::SyntheticSatellite) = get_gnss_system(ss.sat)
@inline get_carrier_to_noise_density_ratio(ss::SyntheticSatellite) =
    get_carrier_to_noise_density_ratio(ss.sat)
