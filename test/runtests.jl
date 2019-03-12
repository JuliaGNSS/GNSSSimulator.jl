using Test, GNSSSimulator, Rotations, CoordinateTransformations, StaticArrays, GNSSSignals, Random, Statistics, Unitful, LinearAlgebra
import Unitful: dB, Hz, kHz, MHz, GHz, dBHz, m, s, rad, °

const STAT_ATT = RotXYZ(0.1, 0.2, 0.3)
const DYN_ATT = [[RotXYZ(0.1, 0.2, 0.3)]; [RotXYZ(0.3, 0.4, 0.5)]; [RotXYZ(0.5, 0.6, 0.7)]]
const STD_ROLL = 0.05
const STD_PITCH = 0.02
const STD_YAW = 0.05
const ΔROLL_PER_S = 0.2rad / 1s
const ΔPITCH_PER_S = 0.2rad / 1s
const ΔYAW_PER_S = 0.2rad / 1s
const CART_COORD = SVector{3}([1.0; 2.0; 3.0])
const SPH_COORD = SphericalFromCartesian()(CART_COORD)
const CART_COORD_VEC = [[CART_COORD]; [2 * CART_COORD]; [3 * CART_COORD]]
const SPH_COORD_VEC = SphericalFromCartesian().(CART_COORD_VEC)
const DYN_DOA_CART = GNSSSimulator.DynamicDOA(CART_COORD_VEC, 1.0Hz)
const DYN_DOA_SPH = GNSSSimulator.DynamicDOA(SPH_COORD_VEC, 1.0Hz)
const LOTHARS_DOAS = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
                     -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
                      0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]
const NUM_ANTS = 4
const EARTH_RADIUS = 6_360_000m

Random.seed!(1234)

include("attitude.jl")
include("doa.jl")
include("existence.jl")
include("gen_example_sat_channels.jl")
include("jammer.jl")
include("measurement.jl")
include("noise.jl")
include("pseudo_post_corr_signal.jl")
include("satellite.jl")
include("system.jl")
