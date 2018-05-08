module GNSSSimulator

  using CoordinateTransformations, StaticArrays, Rotations, PhasedArray
  import Base.transpose

  export init_measurement

  include("general.jl")
  include("phased_array.jl")
end
