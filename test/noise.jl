@testset "Noise" begin

    noise = Noise(randn(ComplexF64), 3.0)
    @inferred propagate(noise, 1μs)
    @inferred get_noise(noise)
    noises = [get_noise(propagate(noise, 1μs)) for i = 1:1000]
    @test noises'noises / 1000 ≈ 9.0 atol = 0.1
end
