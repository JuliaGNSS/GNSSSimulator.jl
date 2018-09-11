@testset "Satellite" begin
    cn0 = 10.0dBHz
    sample_freq = 4e6Hz
    interm_freq = 100_000Hz
    center_freq = 1_575_420_000Hz
    code_length = 1023
    test_range = 1:100
    gnss_system = GPSL1()
    attitude = RotXYZ(0.0, 0.0, 0.0)
    sat = GNSSSimulator.Satellite(prn = 1, enu_doa = CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), CN0 = cn0)
    measurement = @inferred GNSSSimulator.init_sim_emitter_signal(sat, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(100)
    doppler = calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, center_freq)
    code_freq = 1_023_000Hz + doppler * 1_023_000Hz / center_freq
    sat_user_distance = calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    code_phase = GNSSSimulator.calc_code_phase(sat_user_distance, code_freq, code_length)
    carrier_phase = calc_carrier_phase(sat_user_distance, center_freq + doppler)

    @test signal ≈ cis.(2π * (interm_freq + doppler) / sample_freq * test_range + carrier_phase) .* gen_code(gnss_system, test_range, code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0) / sample_freq)

    @test @inferred(GNSSSimulator.calc_code_phase(0m, 1_023_000Hz, 1023)) == 0
    @test @inferred(GNSSSimulator.calc_code_phase(293.255132m, 1_023_000Hz, 1023)) ≈ 1 rtol = 1e-3 # 1 Chip is around 300m 

    @test @inferred(GNSSSimulator.calc_amplitude_from_cn0(45dBHz, 4e6Hz)) ≈ sqrt(10^(45 / 10) / 4e6) 
    @test @inferred(GNSSSimulator.calc_amplitude_from_cn0(42dBHz, 8e6Hz)) ≈ sqrt(10^(42 / 10) / 8e6)
end