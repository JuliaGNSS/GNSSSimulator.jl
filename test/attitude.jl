@testset "Attitude" begin
    @testset "Static Attitude" begin
        @test propagate(RotXYZ(1.0, 2.0, 3.0), 1ms) == RotXYZ(1.0, 2.0, 3.0)

        noisy_static_attitude = NoisyStaticAttitude(RotXYZ(1.0, 2.0, 3.0), 0.1, 0.2, 0.3)
        next_noisy_static_attitude = propagate(noisy_static_attitude, 1ms)
        @test next_noisy_static_attitude.base_attitude == RotXYZ(1.0, 2.0, 3.0)

        @test std([get_attitude(propagate(noisy_static_attitude, 1ms)).theta1 for i = 1:1000]) ≈ 0.1 atol = 0.03
        @test std([get_attitude(propagate(noisy_static_attitude, 1ms)).theta2 for i = 1:1000]) ≈ 0.2 atol = 0.03
        @test std([get_attitude(propagate(noisy_static_attitude, 1ms)).theta3 for i = 1:1000]) ≈ 0.3 atol = 0.03

    end

    @testset "Dynamic Attitude" begin
        attitudes = [RotXYZ(0.1*i, 0.1*i+0.1, 0.1*i+0.2) for i = 0:10]
        dynamic_attitude = DynamicAttitude(attitudes, 0.0s, 1Hz)
        next_dynamic_attitude = propagate(dynamic_attitude, 1s)
        forward_dynamic_attitude = propagate(dynamic_attitude, 20s)
        linear_dynamic_attitude = LinearDynamicAttitude(RotXYZ(0.0, 0.1, 0.2), 0.1rad/s, 0.1rad/s, 0.1rad/s)
        next_linear_dynamic_attitude = propagate(linear_dynamic_attitude, 1s)

        @test next_dynamic_attitude == DynamicAttitude(attitudes, 1.0s, 1Hz)
        @test get_attitude(next_dynamic_attitude) ≈ RotXYZ(0.1, 0.2, 0.3)
        @test get_attitude(forward_dynamic_attitude) == RotXYZ(1.0, 1.1, 1.2)
        @test next_linear_dynamic_attitude.attitude ≈ RotXYZ(0.1, 0.2, 0.3)
        @test get_attitude(next_linear_dynamic_attitude) ≈ RotXYZ(0.1, 0.2, 0.3)

        attitudes2 = [RotXYZ(0.0, 0.1, 0.2) for i = 0:10]
        noisy_dynamic_attitude = NoisyDynamicAttitude(attitudes2, 0.0s, 1Hz, 0.1, 0.2, 0.3)
        next_noisy_dynamic_attitude = propagate(noisy_dynamic_attitude, 1s)
        @test next_noisy_dynamic_attitude.attitudes == attitudes2
        noisy_linear_dynamic_attitude = NoisyLinearDynamicAttitude(RotXYZ(0.0, 0.1, 0.2), 0.0rad/s, 0.0rad/s, 0.0rad/s, 0.1, 0.2, 0.3)
        next_noisy_linear_dynamic_attitude = propagate(noisy_linear_dynamic_attitude, 1s)
        @test next_noisy_linear_dynamic_attitude.base_attitude ≈ RotXYZ(0.0, 0.1, 0.2)

        @test std([get_attitude(propagate(noisy_dynamic_attitude, 1ms)).theta1 for i = 1:1000]) ≈ 0.1 atol = 0.03
        @test std([get_attitude(propagate(noisy_dynamic_attitude, 1ms)).theta2 for i = 1:1000]) ≈ 0.2 atol = 0.03
        @test std([get_attitude(propagate(noisy_dynamic_attitude, 1ms)).theta3 for i = 1:1000]) ≈ 0.3 atol = 0.03

        @test std([get_attitude(propagate(noisy_linear_dynamic_attitude, 1ms)).theta1 for i = 1:1000]) ≈ 0.1 atol = 0.03
        @test std([get_attitude(propagate(noisy_linear_dynamic_attitude, 1ms)).theta2 for i = 1:1000]) ≈ 0.2 atol = 0.03
        @test std([get_attitude(propagate(noisy_linear_dynamic_attitude, 1ms)).theta3 for i = 1:1000]) ≈ 0.3 atol = 0.03
    end
end
