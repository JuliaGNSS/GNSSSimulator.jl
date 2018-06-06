module GNSSSimulator

  using DocStringExtensions, Rotations, JuliennedArrays, Unitful
  import Base.transpose

  export sim_post_corr_measurement, sim_doas, sim_existing_sats, sim_pseudo_post_corr_signal, sim_attitude, sim_noise, sim_gain_phase_mism_and_crosstalk, sim_steering_vectors

  include("general.jl")
  include("phased_array.jl")
end
