@testset "Gain and phase mismatch and crosstalk" begin

    num_ants = 4
    num_sats = 100
    gen_gen_noise = GNSSSimulator.init_gen_gain_and_phase_mism_and_crosstalk(
            num_ants,
            init_phase_mism_betw_ant_std,
            phase_mism_over_time_std,
            init_gain_mism_betw_ant_std,
            gain_mism_over_time_std,
            init_crosstalk_to_direct_power_dB,
            init_crosstalk_ampl_std,
            init_crosstalk_phase_std,
            crosstalk_ampl_over_time_std,
            crosstalk_phase_over_time_std
        )
        
    noise = gen_gen_noise(0, trues(num_sats))
    @test noise * noise' / num_sats ≈ eye(num_ants) atol = 0.4

end
