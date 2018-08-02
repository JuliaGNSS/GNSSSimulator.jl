module GNSSSimulator

  using DocStringExtensions, Rotations, JuliennedArrays, Unitful, CoordinateTransformations, Parameters, GNSSSignals
  import FunctionWrappers
  import Base.transpose
  import FunctionWrappers: FunctionWrapper

  export sim_post_corr_measurement,
    sim_pre_corr_measurement,
    sim_doas,
    sim_existing_sats,
    sim_pseudo_post_corr_signal,
    sim_attitude,
    sim_noise,
    sim_gain_phase_mism_and_crosstalk,
    sim_steering_vectors,
    sim_existing_interfs,
    sim_interf_doas,
    sim_pseudo_post_corr_interf_signal,
    GNSSSystem,
    LinearDynamicDOA,
    Satellite,
    CWJammer,
    LinearDynamicAttitude,
    gpsl1_system,
    gpsl5_system


  include("general.jl")
  include("phased_array.jl")
  include("system.jl")
  include("jammer.jl")
  include("measurement.jl")
  include("satellite.jl")

end
