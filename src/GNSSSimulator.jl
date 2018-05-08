module GNSSSimulator

  using CoordinateTransformations, StaticArrays, Rotations
  import Base.transpose

  export init_measurement

  include("general.jl")
  include("phased_array.jl")
end
