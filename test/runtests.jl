using Base.Test, GNSSSimulator, GNSSSignals, Rotations, CoordinateTransformations, JuliennedArrays

#include("general.jl")
#include("phased_array.jl")
include("jammer.jl")
include("satellite.jl")
include("measurement.jl")