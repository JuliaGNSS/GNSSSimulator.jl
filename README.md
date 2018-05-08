[![pipeline status](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
# GNSSSimulator
Simulate GNSS signals. Currently it only provides pseudo post correlation signals.

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
Pkg.clone("git@git.rwth-aachen.de:nav/GNSSSimulator.jl.git")
```

## Usage

```julia
using GNSSSimulator
measurement = init_measurement(init_measurement(
        attitudes,
        existing_sats,
        sat_doa_carts,
        get_steer_vec)
```

## Todo

This is still missing:

* DOA calculation based on absolute time and NMEA
* PRN simulation (currently it only provides a random post correlation signal for each satellite)
* Multipath
* Jammer
* Spoofer

## Nice to have

* Attitude and trajectory simulation based on Google API
* Multipath simulation based on Google API

## License

MIT License
