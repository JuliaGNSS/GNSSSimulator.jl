const EARTH_RADIUS = 6_360_000m
const SPEED_OF_LIGHT = 299_792_458m / 1s

"""
$(SIGNATURES)

Approximates the distance between user and satellite based on the distance from earth center to satellite `distance_from_earth_center`
and the signal direction of arrival `enu_doa`.
"""
function calc_init_sat_user_distance(distance_from_earth_center, enu_doa)
    init_doa = SphericalFromCartesian()(sim_doa(0s, enu_doa))
    EARTH_RADIUS * cos(init_doa.ϕ + π / 2) + sqrt((EARTH_RADIUS * cos(init_doa.ϕ + π / 2))^2 - EARTH_RADIUS^2 + distance_from_earth_center^2)
end

"""
$(SIGNATURES)

Calculate carrier phase based on the distance between user and satellite `sat_user_distance` and frequency `freq`.
"""
function calc_carrier_phase(sat_user_distance, freq)
    mod2pi(2π * freq * sat_user_distance / SPEED_OF_LIGHT)
end

"""
$(SIGNATURES)

Calculates the doppler based on elevation. Assumes satellite velocity only in elevation direction. Currently there is only a positive doppler.
"""
function calc_init_doppler(distance_from_earth_center, enu_doa, velocity, center_freq)
    init_doa = SphericalFromCartesian()(sim_doa(0s, enu_doa))
    center_freq / SPEED_OF_LIGHT * velocity * cos(π / 2 - asin(EARTH_RADIUS * sin(π / 2 + init_doa.ϕ) / distance_from_earth_center))
end

"""
$(SIGNATURES)

Calculate doppler based on the relative velocity `relative_velocity` and the frequency `center_freq`.
"""
function doppler(relative_velocity, center_freq)
    center_freq / SPEED_OF_LIGHT * relative_velocity
end

"""
$(SIGNATURES)

Create a random complex noise signal.
"""
function gen_noise(num_ants, num_samples)
    complex.(randn(num_samples, num_ants), randn(num_samples, num_ants)) / sqrt(2)
end