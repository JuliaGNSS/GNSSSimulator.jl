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
using GNSSSimulator, Rotations
using GNSSSimulator: Hz, dBHz, GPSL1
gpsl1 = GPSL1()
sample_freq = 2e6Hz
sat = ConstantDopplerSatellite(1, gpsl1, carrier_doppler = 1000.0Hz, carrier_phase = π / 2, code_phase = 100.0, cn0 = 45dBHz)
emitters = (sat,)
receiver = Receiver(1.0, RotXYZ(0.0, 0.0, 0.0), sqrt(sample_freq / 1Hz))
received_signal = ReceivedSignal(emitters, receiver)
measurement1 = get_measurement(received_signal)
next_received_signal = propagate(received_signal, 1/sample_freq)
measurement2 = get_measurement(next_received_signal)
```

## Todo

* DOA calculation based on absolute time and NMEA

## Nice to have

* Attitude and trajectory simulation based on Google API
* Multipath simulation based on Google API

## License

MIT License
