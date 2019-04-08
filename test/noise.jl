@testset "Noise" begin
    noise_power = 10dB
    noise = reduce(hcat, [sim_noise(noise_power, NUM_ANTS) for i = 1:1000])

    # @test all(isapprox.(var(noise, 2, corrected = false), uconvertp(NoUnits, noise_power), atol = 0.5)) # in Julia v0.6
    @test all(isapprox.(var(noise, corrected = false, dims = 2), uconvertp(NoUnits, noise_power), atol = 0.5))

    noise = gen_noise(4, 100000, 4e6Hz)
    @test noise' * noise / 100000 ≈ 4e6I atol = 1e5
end
