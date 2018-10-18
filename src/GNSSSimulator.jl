module GNSSSimulator

    using
        DocStringExtensions,
        Rotations,
        JuliennedArrays,
        Unitful,
        CoordinateTransformations,
        Parameters,
        StaticArrays,
        GNSSSignals,
        FunctionWrappers,
        LinearAlgebra,
        Destruct
        
    import Base.transpose
    import Unitful: Hz, rad, s, m, dB, °, dBHz

    export
        CWJammer,
        SatelliteChannel,
        Satellite,
        DynamicExistence,
        sim_attitude,
        sim_doa,
        sim_existence,
        gen_example_sat_channels,
        calc_amplitude_from_jnr,
        sim_post_corr_measurement,
        init_sim_measurement,
        sim_noise,
        sim_gain_phase_mism_and_crosstalk,
        sim_pseudo_post_corr_signal,
        calc_amplitude_from_cn0,
        calc_code_phase,
        calc_init_sat_user_distance,
        calc_carrier_phase,
        calc_init_doppler,
        doppler,
        gen_noise

    const LOTHARS_DOAS = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
                         -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
                          0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]

    abstract type AbstractDynamicAttitude end
    abstract type AbstractDynamicDOA end
    abstract type AbstractDynamicExistence end
    abstract type AbstractNoisyPseudoPostCorr end
    abstract type AbstractEmitter end
    abstract type AbstractJammer <: AbstractEmitter end

    struct NoisyStaticAttitude
        attitude::RotXYZ{Float64}
        roll_std::Float64
        pitch_std::Float64
        yaw_std::Float64
    end

    struct DynamicAttitude <: AbstractDynamicAttitude
        attitude::Array{RotXYZ{Float64}, 1}
        sample_freq::typeof(1.0Hz)
    end

    @with_kw struct LinearDynamicAttitude <: AbstractDynamicAttitude
        init_attitude::RotXYZ{Float64}
        Δroll::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
        Δpitch::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
        Δyaw::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    end

    struct NoisyDynamicAttitude <: AbstractDynamicAttitude
        attitude::Array{RotXYZ{Float64}, 1}
        sample_freq::typeof(1.0Hz)
        roll_std::Float64
        pitch_std::Float64
        yaw_std::Float64
    end

    @with_kw struct NoisyLinearDynamicAttitude <: AbstractDynamicAttitude
        init_attitude::RotXYZ{Float64}
        Δroll::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
        Δpitch::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
        Δyaw::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
        roll_std::Float64
        pitch_std::Float64
        yaw_std::Float64
    end

    struct DynamicDOA{T <: Union{Vector{SVector{3, R}}, Vector{Spherical{R}}} where R <: Real} <: AbstractDynamicDOA
        doas::T
        sample_freq::typeof(1.0Hz)
    end

    @with_kw struct LinearDynamicDOA{R <: Real} <: AbstractDynamicDOA
        init_DOA::Spherical{R}
        Δazimuth::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
        Δelevation::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    end

    struct DynamicExistence <: AbstractDynamicExistence
        existence::Vector{Bool}
        sample_freq::typeof(1.0Hz)
    end

    struct SatelliteChannel{
        T <: Union{SVector{3, R}, AbstractDynamicDOA} where R <: Real,
        S <: Union{Complex{Float64}, AbstractNoisyPseudoPostCorr},
        E <: Union{Bool, AbstractDynamicExistence},
        IT <: Union{SVector{3, R}, AbstractDynamicDOA} where R <: Real,
        IS <: Union{Complex{Float64}, AbstractNoisyPseudoPostCorr},
        IE <: Union{Bool, AbstractDynamicExistence}
    }
        channel::Int
        enu_doa::T
        signal::S
        exists::E
        interf_enu_doa::IT
        interf_signal::IS
        interf_exists::IE
    end

    struct SatelliteChannelState
        doa::SVector{3, <:Real}
        signal::Complex{Float64}
        exists::Bool
        interf_doa::SVector{3, <:Real}
        interf_signal::Complex{Float64}
        interf_exists::Bool
    end

    struct InternalStates
        sat_channels::Vector{SatelliteChannelState}
        attitude::RotXYZ
        gain_phase_mism_crosstalk::Matrix{Complex{Float64}}
    end

    struct CWJammer{
        T <: Union{SVector{3, R}, AbstractDynamicDOA} where R<:Real,
        E <: Union{Bool, AbstractDynamicExistence}
    } <: AbstractJammer
        id::Int
        enu_doa::T
        relative_velocity::typeof(1.0m / 1.0s)
        JNR::typeof(1.0dB)
        exists::E
    end

    @with_kw struct Satellite{
        T <: Union{SVector{3, R}, AbstractDynamicDOA} where R<:Real,
        E <: Union{Bool, AbstractDynamicExistence}
    } <: AbstractEmitter
        prn::Int
        enu_doa::T
        velocity::typeof(1.0m / 1.0s) = 14_000.0m / 3.6s
        CN0::typeof(1.0dBHz) = 45.0dBHz
        distance_from_earth_center::typeof(1.0m) = 26_560_000.0m
        exists::E = true
    end

    struct NoisyPseudoPostCorr <: AbstractNoisyPseudoPostCorr
        ampl::Float64
        phase::Float64 # in [rad]. Input in [°] will be converted to rad.
        ampl_std::Float64
        phase_std::Float64 # in [rad]. Input in [°] will be converted to rad.
    end

    struct EmitterInternalStates
        doppler::typeof(1.0Hz)
        carrier_phase::Float64 # in [rad]. Input in [°] will be converted to rad.
        code_phase::Float64 # in [rad]. Input in [°] will be converted to rad.
    end

    """
    $(SIGNATURES)

    Generates a default struct of type `SatelliteChannel` if no interference specified. No interference is assumed in that case.
    """
    function SatelliteChannel(channel, enu_doa, signal, exists)
        interf_enu_doa = SVector{3}(0.0, 0.0, 1.0)
        interf_signal = 0.0 + 0.0im
        interf_exists = false
        SatelliteChannel(channel, enu_doa, signal, exists, interf_enu_doa, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)

    Generates a default struct of type `SatelliteChannel` if DOA of static satellite signal given in spherical coordinates. DOA is converted to Cartesian coordinates.
    """
    function SatelliteChannel(channel, enu_doa::Spherical, signal, exists, interf_enu_doa, interf_signal, interf_exists)
        enu_doa_cart = CartesianFromSpherical()(enu_doa)
        SatelliteChannel(channel, enu_doa_cart, signal, exists, interf_enu_doa, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)

    Generates a default struct of type `SatelliteChannel` if DOA of static interference signal given in spherical coordinates. DOA is converted to Cartesian coordinates.
    """
    function SatelliteChannel(channel, enu_doa, signal, exists, interf_enu_doa::Spherical, interf_signal, interf_exists)
        interf_enu_doa_cart = CartesianFromSpherical()(interf_enu_doa)
        SatelliteChannel(channel, interf_enu_doa, signal, exists, interf_enu_doa_cart, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)

    Generates a default struct of type `SatelliteChannel` if DOA of both static satellite and static interference signal given in spherical coordinates. Both DOAs are converted to Cartesian coordinates.
    """
    function SatelliteChannel(channel, enu_doa::Spherical, signal, exists, interf_enu_doa::Spherical, interf_signal, interf_exists)
        enu_doa_cart = CartesianFromSpherical()(enu_doa)
        interf_enu_doa_cart = CartesianFromSpherical()(interf_enu_doa)
        SatelliteChannel(channel, interf_enu_doa, signal, exists, interf_enu_doa_cart, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)

    Creates a struct of type `NoisyPseudoPostCorr` if the satellite signal power or interference power given in [dB]. Assuming no noise.
    """
    function NoisyPseudoPostCorr(signal_power::Unitful.Gain, phase, ampl_std, phase_std)
        NoisyPseudoPostCorr(sqrt(uconvertrp(NoUnits, signal_power)), phase, ampl_std, phase_std)
    end

    include("system.jl")
    include("attitude.jl")
    include("doa.jl")
    include("existence.jl")
    include("gen_example_sat_channels.jl")
    include("jammer.jl")
    include("measurement.jl")
    include("noise.jl")
    include("phased_array_uncertainties.jl")
    include("pseudo_post_corr_signal.jl")
    include("satellite.jl")

end
