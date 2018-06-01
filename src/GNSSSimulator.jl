module GNSSSimulator

  using CoordinateTransformations, StaticArrays, Rotations, JuliennedArrays
  import Base.transpose

  export init_measurement, TemporalData

  include("general.jl")
  include("phased_array.jl")
end
