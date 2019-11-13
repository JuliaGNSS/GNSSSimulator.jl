[![pipeline status](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
# GNSSSimulator
Simulate GNSS signals, jammer, spoofer, multipath.

## Features

 * GNSS signals
 * Jammers
 * Structural interference (spoofer, multipath)
 * Phased arrays:
  * Gain and phase mismatches
  * Crosstalk effects
  * Steering vectors
  * Attitude

## Getting started

Install:
```julia
julia> ]
pkg> add git@git.rwth-aachen.de:nav/GNSSSimulator.jl.git
```

## Usage

```julia
using GNSSSimulator
using GNSSSimulator: GPSL1, Hz
sample_freq = 2e6Hz
sat = ConstantDopplerSatellite(GPSL1, 1)
emitters = (sat,)
receiver = Receiver(sample_freq)
measurement1, next_receiver1, next_emitters1 = get_measurement(2000, receiver, emitters)
measurement2, next_receiver2, next_emitters2 = get_measurement(2000, next_receiver1, next_emitters1)
```

## Todo

* DOA calculation based on absolute time and NMEA

## License

MIT License
