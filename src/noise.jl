"""
$(SIGNATURES)

Simulates noise for one measurement of all antennas of a satellite channel.
"""
function sim_noise(noise_power, num_ants)
    complex.(randn(num_ants), randn(num_ants)) ./ sqrt(2) .* sqrt(uconvertp(NoUnits, noise_power))
end

"""
$(SIGNATURES)

Create a random complex noise signal.
"""
function gen_noise(num_ants, num_samples, sample_freq)
    complex.(randn(num_samples, num_ants), randn(num_samples, num_ants)) / sqrt(2) * sqrt(sample_freq / Hz)
end
