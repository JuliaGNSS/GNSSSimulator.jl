@testset "Post Corr Measurement" begin
    num_sats = 4
    sat_channels = [GNSSSimulator.SatelliteChannel(i, SVector{3}(LOTHARS_DOAS[:,i]), 1 * cis(pi/2), true, SVector{3}(0.0,0.0,1.0), 0.0 + 0.0im, false) for i = 1:num_sats]
    attitudes = STAT_ATT
    gain_phase_mism_and_crosstalk = sim_gain_phase_mism_and_crosstalk(NUM_ANTS, 0.031) # 0.031 ≈ -15dB
    get_steer_vec(doa) = [1 + 0im, 1 + 0im, 1 + 0im, 1 + 0im]
    noise_power = -Inf * 1dB # no noise

    post_corr_measurement = sim_post_corr_measurement(
        sat_channels,
        attitudes,
        gain_phase_mism_and_crosstalk,
        get_steer_vec,
        noise_power)

    𝐘, internal_states = post_corr_measurement(1s)

    @test size(𝐘) == (NUM_ANTS, num_sats)
    @test 𝐘[:, 2] == gain_phase_mism_and_crosstalk(1s) * [1 + 0im, 1 + 0im, 1 + 0im, 1 + 0im] .* (1 * cis(pi/2)) 
    @test size(internal_states.gain_phase_mism_crosstalk) == (NUM_ANTS, NUM_ANTS)
    @test internal_states.sat_channels[1].doa == SVector{3}([0.6409; -0.6409; 0.4226])
    @test internal_states.sat_channels[1].signal == 1 * cis(pi/2)
    @test internal_states.sat_channels[1].exists == true
    @test internal_states.sat_channels[1].interf_doa == SVector{3}(0.0, 0.0, 1.0)
    @test internal_states.sat_channels[1].interf_signal == 0.0 + 0.0im
    @test internal_states.sat_channels[1].interf_exists == false
    @test internal_states.attitude == RotXYZ(0.1, 0.2, 0.3)
end

@testset "Continuity of measurement" begin 
    sample_freq = 4e6Hz
    interm_freq = 100_000Hz
    gnss_system = GNSSSimulator.GPSL1()
    jammer = GNSSSimulator.CWJammer(1, CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), 0.0m / 1.0s, 20.0dB, true)
    sat = GNSSSimulator.Satellite(prn = 1, enu_doa = CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)))
    attitude = RotXYZ(0.0, 0.0, 0.0)
    emitters = [sat, jammer]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal1, internal_states = measurement(10)
    next_measurement, signal_short1, internal_states = measurement(5)
    next_measurement, signal_short2, internal_states = next_measurement(5)
    signal2 = [signal_short1; signal_short2]
    @test signal1 ≈ signal2
end

@testset "Measurement with one Jammer" begin 
    sample_freq = 4e6Hz
    interm_freq = 100_000Hz
    jammer_power = 20.0dB
    gnss_system = GNSSSimulator.GPSL1()
    jammer = GNSSSimulator.CWJammer(1, CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), 0.0m / 1s, jammer_power, true)
    attitude = RotXYZ(0.0, 0.0, 0.0)
    emitters = [jammer]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement = @inferred GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal, internal_states = measurement(100)
    @test signal ≈ transpose(get_steer_vec(Spherical(0.0, 0.0, 1.0))) .* cis.(2π * interm_freq / sample_freq * (1:100)) * uconvertrp(NoUnits, jammer_power)

    jammer = GNSSSimulator.CWJammer(1, CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), 0.0m / 1s, jammer_power, false)
    measurement = @inferred GNSSSimulator.init_sim_measurement([jammer], gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal, internal_states = measurement(100)
    @test signal ≈ zeros(100, 4)
end