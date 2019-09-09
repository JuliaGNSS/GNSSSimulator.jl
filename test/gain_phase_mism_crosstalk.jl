@testset "Gain and phase mismatch and crosstalk" begin

    gain_and_phase_mism_crosstalk = SMatrix{2,2}(1,2,3,4)
    @test propagate(gain_and_phase_mism_crosstalk, 1ms) == gain_and_phase_mism_crosstalk
end
