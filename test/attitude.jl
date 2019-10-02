@testset "Attitude" begin
    @testset "Static Attitude" begin
        @test GNSSSimulator.propagate(RotXYZ(1.0, 2.0, 3.0), 1ms, Random.GLOBAL_RNG) == RotXYZ(1.0, 2.0, 3.0)

        noisy_static_attitude = @inferred NoisyStaticAttitude(RotXYZ(1.0, 2.0, 3.0), 0.1, 0.2, 0.3)
        rng = MersenneTwister(1234)
        next_noisy_static_attitude = @inferred GNSSSimulator.propagate(noisy_static_attitude, 1ms, rng)
        @test next_noisy_static_attitude.base_attitude == RotXYZ(1.0, 2.0, 3.0)
        rng = MersenneTwister(1234)
        @test get_attitude(next_noisy_static_attitude) == RotXYZ(1.0 + randn(rng) * 0.1, 2.0 + randn(rng) * 0.2, 3.0 + randn(rng) * 0.3)

    end

    @testset "Dynamic Attitude" begin
        attitudes = [RotXYZ(0.1*i, 0.1*i+0.1, 0.1*i+0.2) for i = 0:10]
        dynamic_attitude = @inferred DynamicAttitude(attitudes, 0.0s, 1Hz)
        next_dynamic_attitude = @inferred GNSSSimulator.propagate(dynamic_attitude, 1s, Random.GLOBAL_RNG)
        forward_dynamic_attitude = @inferred GNSSSimulator.propagate(dynamic_attitude, 20s, Random.GLOBAL_RNG)
        linear_dynamic_attitude = @inferred LinearDynamicAttitude(RotXYZ(0.0, 0.1, 0.2), 0.1rad/s, 0.1rad/s, 0.1rad/s)
        next_linear_dynamic_attitude = @inferred GNSSSimulator.propagate(linear_dynamic_attitude, 1s, Random.GLOBAL_RNG)

        @test next_dynamic_attitude == DynamicAttitude(attitudes, 1.0s, 1Hz)
        @test get_attitude(next_dynamic_attitude) ≈ RotXYZ(0.1, 0.2, 0.3)
        @test get_attitude(forward_dynamic_attitude) == RotXYZ(1.0, 1.1, 1.2)
        @test next_linear_dynamic_attitude.attitude ≈ RotXYZ(0.1, 0.2, 0.3)
        @test get_attitude(next_linear_dynamic_attitude) ≈ RotXYZ(0.1, 0.2, 0.3)

        attitudes2 = [RotXYZ(0.0, 0.1, 0.2) for i = 0:10]
        noisy_dynamic_attitude = @inferred NoisyDynamicAttitude(attitudes2, 0.0s, 1Hz, 0.1, 0.2, 0.3)

        rng = MersenneTwister(1234)
        next_noisy_dynamic_attitude = @inferred GNSSSimulator.propagate(noisy_dynamic_attitude, 1s, rng)
        @test next_noisy_dynamic_attitude.attitudes == attitudes2
        rng = MersenneTwister(1234)
        @test get_attitude(next_noisy_dynamic_attitude) == RotXYZ(0.0 + randn(rng) * 0.1, 0.1 + randn(rng) * 0.2, 0.2 + randn(rng) * 0.3)
        noisy_linear_dynamic_attitude = @inferred NoisyLinearDynamicAttitude(RotXYZ(0.0, 0.1, 0.2), 0.1rad/s, 0.2rad/s, 0.3rad/s, 0.1, 0.2, 0.3)
        rng = MersenneTwister(1234)
        next_noisy_linear_dynamic_attitude = @inferred GNSSSimulator.propagate(noisy_linear_dynamic_attitude, 1s, rng)
        @test next_noisy_linear_dynamic_attitude.base_attitude ≈ RotXYZ(0.1, 0.3, 0.5)
        rng = MersenneTwister(1234)
        @test get_attitude(next_noisy_linear_dynamic_attitude) ≈ RotXYZ(0.1 + randn(rng) * 0.1, 0.3 + randn(rng) * 0.2, 0.5 + randn(rng) * 0.3)
    end
end
