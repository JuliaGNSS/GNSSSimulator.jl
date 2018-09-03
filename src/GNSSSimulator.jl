module GNSSSimulator

    using DocStringExtensions, Rotations, JuliennedArrays, Unitful, CoordinateTransformations, Parameters, StaticArrays
    import Base.transpose
    import Unitful: s, Hz, dB, rad, °
        
    export sim_post_corr_measurement,
        sim_doa,
        sim_existence,
        sim_pseudo_post_corr_signal,
        sim_attitude,
        sim_noise,
        sim_attitude,
        sim_gain_phase_mism_and_crosstalk,
        sim_steering_vectors,
        gen_example_sat_channels

    const LOTHARS_DOAS = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
                         -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
                          0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]    

    abstract type AbstractDynamicDOA end
    abstract type AbstractDynamicExistence end

    struct SatelliteChannelState
        doa::SVector{3, <:Real}
        signal::Complex{Float64}
        exists::Bool
        interf_doa::SVector{3, <:Real}
        interf_signal::Complex{Float64}
        interf_exists::Bool
    end
        
    struct DynamicDOA{T <: Union{Vector{SVector{3, R}}, Vector{Spherical{R}}} where R <: Real, F <: Unitful.Frequency} <: AbstractDynamicDOA
        doas::T
        sample_freq::F
    end
        
    @with_kw struct LinearDynamicDOA{R <: Real, FA <: Unitful.Frequency, FE <: Unitful.Frequency} <: AbstractDynamicDOA
        init_DOA::Spherical{R}
        Δazimuth::FA = 0rad / 1s
        Δelevation::FE = 0rad / 1s
    end
        
    struct DynamicExistence{F <: Unitful.Frequency} <: AbstractDynamicExistence
        existence::Vector{Bool}
        sample_freq::F
    end
        
    struct NoisyPseudoPostCorr
        ampl::Float64
        phase::Float64 # in [rad]
        ampl_std::Float64
        phase_std::Float64 # in [rad]
    end
        
    struct StaticAttitudes
        attitude::RotXYZ{Float64}
    end
        
    struct NoisyStaticAttitudes
        attitude::RotXYZ{Float64}
        roll_std::Float64
        pitch_std::Float64
        yaw_std::Float64
    end
        
    struct DynamicAttitudes{F <: Unitful.Frequency}
        attitude::Array{RotXYZ{Float64}, 1}
        sample_freq::F
    end

    @with_kw struct LinearDynamicAttitudes{FR <: Unitful.Frequency, FP <: Unitful.Frequency, FY <: Unitful.Frequency}
        init_attitude::RotXYZ{Float64}
        Δroll::FR = 0rad / 1s
        Δpitch::FP = 0rad / 1s
        Δyaw::FY = 0rad / 1s
    end
        
    struct NoisyDynamicAttitudes{F <: Unitful.Frequency}
        attitude::Array{RotXYZ{Float64}, 1}
        sample_freq::F
        roll_std::Float64
        pitch_std::Float64
        yaw_std::Float64
    end

    @with_kw struct NoisyLinearDynamicAttitudes{FR <: Unitful.Frequency, FP <: Unitful.Frequency, FY <: Unitful.Frequency}
        init_attitude::RotXYZ{Float64}
        Δroll::FR = 0rad / 1s
        Δpitch::FP = 0rad / 1s
        Δyaw::FY = 0rad / 1s
        roll_std::Float64
        pitch_std::Float64
        yaw_std::Float64
    end

    struct SatelliteChannel{
        T <: Union{SVector{3, R}, AbstractDynamicDOA} where R <: Real,
        S <: Union{Complex{Float64}, NoisyPseudoPostCorr},
        E <: Union{Bool, AbstractDynamicExistence},
        IT <: Union{SVector{3, R}, AbstractDynamicDOA} where R <: Real,
        IS <: Union{Complex{Float64}, NoisyPseudoPostCorr},
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

    struct InternalStates
        sat_channels::Vector{SatelliteChannelState}
        attitude::RotXYZ
        gain_phase_mism_crosstalk::Matrix{Complex{Float64}}
    end
        
    """
    $(SIGNATURES)
        
    Generates a default struct of type "Satellite Channel" if no interference specified. "No interference" is assumed in that case.
    """
    function SatelliteChannel(channel, enu_doa, signal, exists)
        interf_enu_doa = SVector{3}(0.0,0.0,1.0)
        interf_signal = 0.0 + 0.0im
        interf_exists = false 
        SatelliteChannel(channel, enu_doa, signal, exists, interf_enu_doa, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)
        
    Generates a default struct of type "Satellite Channel" if DOA of static satellite signal given in spherical coordinates. DOA is converted to Cartesian coordinates.
    """
    function SatelliteChannel(channel, enu_doa::Spherical, signal, exists, interf_enu_doa, interf_signal, interf_exists)
        enu_doa_cart = CartesianFromSpherical()(enu_doa)
        SatelliteChannel(channel, enu_doa_cart, signal, exists, interf_enu_doa, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)
        
    Generates a default struct of type "Satellite Channel" if DOA of static interference signal given in spherical coordinates. DOA is converted to Cartesian coordinates.
    """
    function SatelliteChannel(channel, enu_doa, signal, exists, interf_enu_doa::Spherical, interf_signal, interf_exists)
        interf_enu_doa_cart = CartesianFromSpherical()(interf_enu_doa)
        SatelliteChannel(channel, interf_enu_doa, signal, exists, interf_enu_doa_cart, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)

    Generates a default struct of type "Satellite Channel" if DOA of both static satellite and static interference signal given in spherical coordinates. Both DOAs are converted to Cartesian coordinates.
    """
    function SatelliteChannel(channel, enu_doa::Spherical, signal, exists, interf_enu_doa::Spherical, interf_signal, interf_exists)
        enu_doa_cart = CartesianFromSpherical()(enu_doa)
        interf_enu_doa_cart = CartesianFromSpherical()(interf_enu_doa)
        SatelliteChannel(channel, interf_enu_doa, signal, exists, interf_enu_doa_cart, interf_signal, interf_exists)
    end

    """
    $(SIGNATURES)

    Creates a struct of type "NoisyPseudoPostCorr" if the satellite signal power or interference power given in dB. Assuming no noise.
    """
    function NoisyPseudoPostCorr(signal_power::Unitful.Gain, phase, ampl_std, phase_std)
        NoisyPseudoPostCorr(sqrt(uconvertp(NoUnits, signal_power)), phase, ampl_std, phase_std)
    end

    """
    $(SIGNATURES)

    Generates a struct of type "StaticAttitudes" if initial values for yaw, pitch and roll angles (in [rad]) are given.
    """
    function StaticAttitudes(init_roll, init_pitch, init_yaw)
        StaticAttitudes(RotXYZ(init_roll, init_pitch, init_yaw))
    end

    include("gen_example_sat_channels.jl")
    include("sim_doa.jl")
    include("sim_existence.jl")
    include("sim_pseudo_post_corr_signal.jl")
    include("sim_noise.jl")
    include("sim_attitude.jl")
    include("phased_array_uncertainties.jl")
    include("measurement.jl")
end