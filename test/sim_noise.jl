@testset "Noise" begin
    num_ants = 4
    noise_power = 10dB
    noise = reduce(hcat, [sim_noise(noise_power, num_ants) for i = 1:1000])

    @test all(isapprox.(var(noise, 2, corrected = false), uconvertp(NoUnits, noise_power), atol = 0.5))
    # @test all(var(noise, 2, corrected = false) .≈ uconvertp(NoUnits, noise_power) atol = 0.5)
end