@testset "Gain and phase mismatch and crosstalk" begin

    @testset "Crosstalk amplitude"
        crosstalk_dB = -10
        gp_mism_crosstalk = GNSSSimulator.init_gen_gain_and_phase_mism_and_crosstalk(2, 0, 0, 0, 0, crosstalk_dB, 0, 0, 0, 0)

        @test abs(gp_mism_crosstalk[1,1]) == abs(gp_mism_crosstalk[2,2])
        @test abs(gp_mism_crosstalk[1,1]) ≈ abs(gp_mism_crosstalk[1,2]) * 10^(crosstalk_dB / 10)
        @test abs(gp_mism_crosstalk[1,2]) == abs(gp_mism_crosstalk[2,1])

    end

end
