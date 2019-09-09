using Test, GNSSSimulator, Rotations, CoordinateTransformations, StaticArrays, GNSSSignals, Random, Statistics, Unitful, LinearAlgebra
import Unitful: dB, Hz, kHz, MHz, GHz, dBHz, m, s, rad, °, ms, μs


const NUM_ANTS = 4
const EARTH_RADIUS = 6_360_000m

Random.seed!(1234)

include("attitude.jl")
include("doa.jl")
include("existence.jl")
include("satellite.jl")
include("jammer.jl")
include("structural_interference.jl")
include("gain_phase_mism_crosstalk.jl")
include("receiver.jl")
include("measurement.jl")
