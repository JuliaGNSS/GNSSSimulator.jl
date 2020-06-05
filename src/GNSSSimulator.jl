module GNSSSimulator

    using
        DocStringExtensions,
        Rotations,
        Unitful,
        CoordinateTransformations,
        StaticArrays,
        GNSSSignals,
        LinearAlgebra,
        Parameters,
        PhasedArray,
        Random,
        StructArrays

    import Base.transpose
    import Unitful: Hz, rad, s, m, dB, °, dBHz, upreferred
    import PhasedArray.get_steer_vec

    const EARTH_RADIUS = 6_360_000m
    const SPEED_OF_LIGHT = 299_792_458m/s

export
    NoisyStaticAttitude,
    DynamicAttitude,
    NoisyDynamicAttitude,
    LinearDynamicAttitude,
    NoisyLinearDynamicAttitude,
    DynamicDOA,
    LinearDynamicDOA,
    DynamicExistence,
    ConstantDopplerSatellite,
    ConstantDopplerStructuralInterference,
    SyntheticSatellite,
    AsymptoticGainPhaseMismCrosstalk,
    CWJammer,
    NoiseJammer,
    Noise,
    Receiver,
    get_attitude,
    get_doa,
    get_existence,
    get_gain_phase_mism_crosstalk,
    get_amplitude,
    get_carrier_doppler,
    get_carrier_phase,
    get_code_doppler,
    get_code_phase,
    get_gnss_system,
    get_carrier_to_noise_density_ratio,
    get_noise_density,
    get_noise_std,
    get_prn,
    get_id,
    get_jammer_to_noise_ratio,
    get_sample_frequency,
    get_intermediate_frequency,
    get_measurement!,
    get_measurement

    propagate(val::Number, Δt, rng) = val
    const cart2sph = SphericalFromCartesian()
    const sph2cart = CartesianFromSpherical()

    function get_sampled_value(time, sample_freq, values)
        index = floor(Int, time * sample_freq) + 1
        index < size(values, 1) ? values[index] : values[end]
    end

    include("attitude.jl")
    include("doa.jl")
    include("existence.jl")
    include("emitter.jl")
    include("satellite.jl")
    include("structural_interference.jl")
    include("synthetic_satellite.jl")
    include("jammer.jl")
    include("gain_phase_mism_crosstalk.jl")
    include("receiver.jl")
    include("measurement.jl")
end
