@testset "Satellite L1" begin
    cn0 = 10.0dBHz
    sample_freq = 4e6Hz
    interm_freq = 100_000Hz
    center_freq = 1_575_420_000Hz
    code_length = 1023
    num_samples = 4000
    test_range = 1:num_samples
    gnss_system = GPSL1()
    attitude = RotXYZ(0.0, 0.0, 0.0)
    sat = Satellite(prn = 1, enu_doa = CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), CN0 = cn0)
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement, init_internal_states = GNSSSimulator.init_sim_measurement([sat], gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal, internal_states = @inferred measurement(num_samples)
    doppler = calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, center_freq)
    code_freq = 1_023_000Hz + doppler * 1_023_000Hz / center_freq
    sat_user_distance = calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    code_phase = GNSSSimulator.calc_code_phase(sat_user_distance, code_freq, code_length)
    carrier_phase = GNSSSimulator.calc_carrier_phase(sat_user_distance, center_freq + doppler)

    @test init_internal_states[1].doppler ≈ doppler
    @test init_internal_states[1].carrier_phase ≈ carrier_phase
    @test init_internal_states[1].code_phase ≈ code_phase

    @test signal ≈ [0.5, 0.5, 0.5, 0.5]' .* cis.(2π * (interm_freq + doppler) / sample_freq * test_range .+ carrier_phase) .* gen_code.(Ref(gnss_system), test_range, code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0) * 1/1Hz)

    next_measurement, signal, internal_states = @inferred next_measurement(num_samples)
    @test signal ≈ [0.5, 0.5, 0.5, 0.5]' .* cis.(2π .* (interm_freq .+ doppler) ./ sample_freq .* (num_samples + 1:2 * num_samples) .+ carrier_phase) .* gen_code.(Ref(gnss_system), (num_samples + 1:2 * num_samples), code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0) * 1/1Hz)
end

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
    attitude = GNSSSimulator.RotXYZ(0.0, 0.0, 0.0)
    emitters = [sat, jammer]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement, init_internal_states = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

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
    measurement, init_internal_states = @inferred GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal, internal_states = measurement(100)
    @test signal ≈ transpose(get_steer_vec(Spherical(0.0, 0.0, 1.0))) .* cis.(2π * interm_freq / sample_freq * (1:100)) * sqrt(uconvertp(NoUnits, jammer.JNR) * sample_freq / 1Hz)
    jammer = GNSSSimulator.CWJammer(1, CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), 0.0m / 1s, jammer_power, false)
    measurement, init_internal_states = @inferred GNSSSimulator.init_sim_measurement([jammer], gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal, internal_states = measurement(100)
    @test signal ≈ zeros(100, 4)
end
