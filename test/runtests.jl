# "RUN": include("GNSSSimulator.jl/test/runtests.jl")
using Base.Test, GNSSSimulator, Rotations, CoordinateTransformations, JuliennedArrays, StaticArrays
using Unitful: dB, Hz, s, rad, °

srand(1234)

include("sim_doa.jl")
include("sim_existence.jl")
include("sim_pseudo_post_corr_signal.jl")
include("sim_noise.jl")
include("sim_attitude.jl")
include("gen_example_sat_channels.jl")
include("GNSSSimulator.jl")
include("measurement.jl")
