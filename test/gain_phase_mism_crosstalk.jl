@testset "Gain and phase mismatch and crosstalk" begin

    @testset "Static" begin
        gain_and_phase_mism_crosstalk = SMatrix{2,2}(1,2,3,4)
        @test propagate(gain_and_phase_mism_crosstalk, 1ms, Random.GLOBAL_RNG) == gain_and_phase_mism_crosstalk

        @test get_gain_phase_mism_crosstalk(gain_and_phase_mism_crosstalk) == SMatrix{2,2}(1,2,3,4)
    end

    @testset "Asymptotic" begin
        a = SVector(1/2, -1/2)
        b = SVector(3/4, 5/4)
        gain_phase_mism_crosstalk = @inferred AsymptoticGainPhaseMismCrosstalk(a, b, e = 2)
        @test gain_phase_mism_crosstalk.a == SVector(1/2, -1/2)
        @test gain_phase_mism_crosstalk.b == SVector(3/4, 5/4)
        @test gain_phase_mism_crosstalk.e == 2
        @test gain_phase_mism_crosstalk.C == SMatrix{2,2,ComplexF64}(I)
        @test gain_phase_mism_crosstalk.t == 0.0s

        gpmc = @inferred get_gain_phase_mism_crosstalk(gain_phase_mism_crosstalk)
        @test gpmc == Diagonal(cis.([3/4; 5/4]))

        next_gain_phase_mism_crosstalk = @inferred GNSSSimulator.propagate(gain_phase_mism_crosstalk, 1000000s, Random.GLOBAL_RNG)
        next_gpmc = @inferred get_gain_phase_mism_crosstalk(next_gain_phase_mism_crosstalk)
        @test next_gpmc ≈ Diagonal(cis.([1/2, -1/2]))
    end
end
