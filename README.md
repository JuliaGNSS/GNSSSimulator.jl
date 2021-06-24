[![Tests](https://github.com/JuliaGNSS/GNSSSimulator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaGNSS/GNSSSimulator.jl/actions)
[![codecov](https://codecov.io/gh/JuliaGNSS/GNSSSimulator.jl/branch/master/graph/badge.svg?token=KCFJHJ4Q2T)](https://codecov.io/gh/JuliaGNSS/GNSSSimulator.jl)
# GNSSSimulator
Simulate GNSS signals, jammers and structural interference.

## Features

  * GNSS signals
  * Jammers
  * Structural interference
  * Phased arrays:
    * Gain and phase mismatches
    * Crosstalk effects
    * Steering vectors
    * Attitude

## Getting started

Install:
```julia
julia> ]
pkg> add GNSSSimulator
```

## Usage

```julia
using GNSSSimulator
using GNSSSimulator: GPSL1, Hz
sample_freq = 2e6Hz
sat = ConstantDopplerSatellite(GPSL1(), 1)
emitters = (sat,)
receiver = Receiver(sample_freq)
measurement1, next_receiver1, next_emitters1 = get_measurement(2000, receiver, emitters)
measurement2, next_receiver2, next_emitters2 = get_measurement(2000, next_receiver1, next_emitters1)
```

## License

MIT License
