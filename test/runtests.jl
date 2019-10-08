using Test, GNSSSimulator, Rotations, CoordinateTransformations, StaticArrays, GNSSSignals, Random, Unitful, LinearAlgebra, PhasedArray
import Unitful: dB, Hz, kHz, MHz, GHz, dBHz, m, s, rad, °, ms, μs

const EARTH_RADIUS = 6_360_000m

include("attitude.jl")
include("doa.jl")
include("existence.jl")
include("satellite.jl")
include("jammer.jl")
include("structural_interference.jl")
include("gain_phase_mism_crosstalk.jl")
include("receiver.jl")
include("measurement.jl")
