@testset "DOAs" begin
    @testset "Static DOAs" begin
        @test propagate(SVector(1.0, 2.0, 3.0), 1s) == SVector(1.0, 2.0, 3.0)
        @test get_doa(SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0, 3.0)

        @test propagate(Spherical(1.0, 0.0, π / 2), 1s) == Spherical(1.0, 0.0, π / 2)
        @test get_doa(Spherical(1.0, 0.0, π / 2)) ≈ SVector(0.0, 0.0, 1.0)
    end

    @testset "Dynamic DOAs" begin
        doas = [SVector(0.0, 1.0+0.1i, 2.0+0.1i) for i = 0:10]
        dynamic_doa = @inferred DynamicDOA(doas, 0.0s, 1.0Hz)
        next_dynamic_doa = @inferred propagate(dynamic_doa, 1s)
        @test next_dynamic_doa == DynamicDOA(doas, 1.0s, 1.0Hz)
        @test get_doa(next_dynamic_doa) == SVector(0.0, 1.1, 2.1)
        @test get_doa(propagate(dynamic_doa, 100s)) ≈ SVector(0.0, 2.0, 3.0)

        linear_dynamic_doa = @inferred LinearDynamicDOA(Spherical(1.0, 0.0, 0.0), 0.1rad/s, 0.0rad/s)
        next_linear_dynamic_doa = @inferred propagate(linear_dynamic_doa, 1s)
        @test get_doa(next_linear_dynamic_doa) ≈ Spherical(1.0, 0.1, 0.0)
    end
end
