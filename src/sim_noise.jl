"""
$(SIGNATURES)

Simulates noise for one measurement (i.e., for all antennas) of a satellite channel.
"""
function sim_noise(noise_power, num_ants)
    complex.(randn(num_ants), randn(num_ants)) ./ sqrt(2) .* sqrt(uconvertp(NoUnits, noise_power))
end