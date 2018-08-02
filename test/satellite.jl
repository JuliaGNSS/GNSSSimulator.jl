@testset "Satellite" begin
    cn0_dB = 10.0
    sample_freq = 4e6
    interm_freq = 100_000
    f₀ = 1_575_420_000
    code_length = 1023
    test_range = 1:100
    gnss_system = GNSSSimulator.gpsl1_system()
    calc_sampled_code, calc_next_code_phase = init_gpsl1_codes()
    attitude = RotXYZ(0.0, 0.0, 0.0)
    sat = Satellite(svid = 1, enu_doa = Spherical(0.0, 0.0, 1.0), CN0 = cn0_dB)
    measurement = @inferred GNSSSimulator.init_sim_emitter_signal(sat, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(100)
    doppler = GNSSSimulator.calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, f₀)
    code_freq = 1_023_000 + doppler * 1_023_000 / f₀
    sat_user_distance = GNSSSimulator.calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    code_phase = GNSSSimulator.calc_code_phase(sat_user_distance, code_freq, code_length)
    carrier_phase = GNSSSimulator.calc_carrier_phase(sat_user_distance, f₀ + doppler)
    @test signal ≈ cis.(2π * (interm_freq + doppler) / sample_freq * test_range + carrier_phase) .* calc_sampled_code(test_range, code_freq, code_phase, sample_freq, 1) .* sqrt(10^(cn0_dB / 10) / sample_freq)

    @test @inferred(GNSSSimulator.calc_code_phase(0, 1_023_000, 1023)) == 0
    @test @inferred(GNSSSimulator.calc_code_phase(293.255132, 1_023_000, 1023)) ≈ 1 rtol = 1e-3 # 1 Chip is around 300m

    @test @inferred(GNSSSimulator.calc_amplitude_from_cn0(45, 1)) ≈ 10^(45 / 20)
    @test @inferred(GNSSSimulator.calc_amplitude_from_cn0(45, 2)) ≈ 10^(42 / 20) atol = 0.2
end