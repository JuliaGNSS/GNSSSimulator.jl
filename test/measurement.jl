@testset "Measurement" begin
    num_sats = 4
    sat_channels = [GNSSSimulator.SatelliteChannel(i, SVector{3}(LOTHARS_DOAS[:,i]), 1 * cis(pi/2), true, SVector{3}(0.0,0.0,1.0), 0.0 + 0.0im, false) for i = 1:num_sats]
    attitudes = GNSSSimulator.StaticAttitudes(RotXYZ(0.1, 0.2, 0.3))
    gain_phase_mism_and_crosstalk = sim_gain_phase_mism_and_crosstalk(NUM_ANTS, 0.031)
    get_steer_vec(doa) = [1 + 0im, 1 + 0im, 1 + 0im, 1 + 0im]
    noise_power = -Inf * 1dB # no noise

    measurement = sim_post_corr_measurement(
        sat_channels,
        attitudes,
        gain_phase_mism_and_crosstalk,
        get_steer_vec,
        noise_power)

    𝐘, internal_states = measurement(1s)

    @test size(𝐘) == (NUM_ANTS, num_sats)
    @test 𝐘[:, 2] == gain_phase_mism_and_crosstalk(1s) * [1 + 0im, 1 + 0im, 1 + 0im, 1 + 0im] .* (1 * cis(pi/2)) 
    @test size(internal_states.gain_phase_mism_crosstalk) == (NUM_ANTS, NUM_ANTS)
    @test internal_states.sat_channels[1].doa == SVector{3}([0.6409; -0.6409; 0.4226])
    @test internal_states.sat_channels[1].signal == 1 * cis(pi/2)
    @test internal_states.sat_channels[1].exists == true
    @test internal_states.sat_channels[1].interf_doa == SVector{3}(0.0,0.0,1.0)
    @test internal_states.sat_channels[1].interf_signal == 0.0 + 0.0im
    @test internal_states.sat_channels[1].interf_exists == false
    @test internal_states.attitude == RotXYZ(0.1, 0.2, 0.3)
end