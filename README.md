[![pipeline status](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/GNSSSimulator.jl/commits/master)
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
doas_over_time = sim_doas()
existing_sats_over_time = sim_existing_sats(trues(11))
pseudo_post_corr_signal_over_time = sim_pseudo_post_corr_signal(11, 0)
attitude_over_time = sim_attitude(0.0, 0.0, 0.0)
gain_phase_mism_and_crosstalk_over_time = sim_gain_phase_mism_and_crosstalk(4, -15)
steering_vectors_over_time = sim_steering_vectors(a -> [a[1] + 0.0im, a[1] + 0.0im, a[2] + 0.0im, a[3] + 0.0im])
noise_over_time = sim_noise(-15, 4)

measurement = sim_pseudo_post_corr_measurement(
    existing_sats_over_time,
    pseudo_post_corr_signal_over_time,
    attitude_over_time,
    doas_over_time,
    gain_phase_mism_and_crosstalk_over_time,
    steering_vectors_over_time,
    noise_over_time)

𝐘, internal_states = measurement(0)
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
