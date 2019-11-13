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
import GNSSSimulator: GPSL1
import Unitful: ms
sample_freq = 2e6Hz
sat = ConstantDopplerSatellite(GPSL1, 1)
emitters = (sat,)
receiver = Receiver(sample_freq)
received_signal = ReceivedSignal(emitters, receiver)
measurement1 = get_measurement(2000, received_signal)
next_received_signal = propagate(received_signal, 1ms)
measurement2 = get_measurement(2000, next_received_signal)
```

## Todo

* DOA calculation based on absolute time and NMEA

## License

MIT License
