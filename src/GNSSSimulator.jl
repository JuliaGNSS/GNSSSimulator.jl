module GNSSSimulator

    using
        DocStringExtensions,
        Rotations,
        Unitful,
        CoordinateTransformations,
        StaticArrays,
        GNSSSignals,
        LinearAlgebra,
        Parameters

    import Base.transpose
    import Unitful: Hz, rad, s, m, dB, °, dBHz

    const EARTH_RADIUS = 6_360_000m
    const SPEED_OF_LIGHT = 299_792_458m/s

export
    propagate,
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
    CWJammer,
    NoiseJammer,
    Receiver,
    ReceivedSignal,
    get_attitude,
    get_doa,
    get_existence,
    get_signal,
    get_gain_phase_mism_crosstalk,
    get_measurement

    const LOTHARS_DOAS = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
                         -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
                          0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]

    propagate(val::Number, Δt) = val
    cart2sph = SphericalFromCartesian()

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
    include("jammer.jl")
    include("gain_phase_mism_crosstalk.jl")
    include("receiver.jl")
    include("measurement.jl")
end
