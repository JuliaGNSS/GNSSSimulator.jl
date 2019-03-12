[![pipeline status](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
# GNSSSimulator
Simulate GNSS signals, jammer, spoofer, multipath.

## Features

 * Phased array simulation:
  * Gain and phase mismatches
  * Crosstalk effects
  * Steering vectors
  * Attitude
 * Pseudo satellite amplitude and phase
 * Standard deviation for variation over time

## Getting started

Install:
```julia
julia> ]
pkg> add git@git.rwth-aachen.de:nav/GNSSSimulator.jl.git
```

## Usage

```julia
using GNSSSimulator, CoordinateTransformations
using GNSSSimulator: Hz, m, s, dB
sample_freq = 4e6Hz
interm_freq = 100_000Hz
gnss_system = GNSSSimulator.GPSL1()
jammer = GNSSSimulator.CWJammer(1, CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), 0.0m / 1.0s, 20.0dB, true)
sat = GNSSSimulator.Satellite(prn = 1, enu_doa = CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)))
attitude = GNSSSimulator.RotXYZ(0.0, 0.0, 0.0)
emitters = [sat, jammer]
get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
measurement, init_internal_states = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

next_measurement, signal1, internal_states = measurement(4000)
```

## Todo

* DOA calculation based on absolute time and NMEA

## Nice to have

* Attitude and trajectory simulation based on Google API
* Multipath simulation based on Google API

## License

MIT License
