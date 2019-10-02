@testset "Gain and phase mismatch and crosstalk" begin

    gain_and_phase_mism_crosstalk = SMatrix{2,2}(1,2,3,4)
    @test propagate(gain_and_phase_mism_crosstalk, 1ms, Random.GLOBAL_RNG) == gain_and_phase_mism_crosstalk

    @test get_gain_phase_mism_crosstalk(gain_and_phase_mism_crosstalk) == SMatrix{2,2}(1,2,3,4)
end
