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
    measurement, init_internal_states = @inferred GNSSSimulator.init_sim_emitter_signal(sat, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(num_samples)
    doppler = calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, center_freq)
    code_freq = 1_023_000Hz + doppler * 1_023_000Hz / center_freq
    sat_user_distance = calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    code_phase = GNSSSimulator.calc_code_phase(sat_user_distance, code_freq, code_length)
    carrier_phase = GNSSSimulator.calc_carrier_phase(sat_user_distance, center_freq + doppler)

    @test init_internal_states.doppler ≈ doppler
    @test init_internal_states.carrier_phase ≈ carrier_phase
    @test init_internal_states.code_phase ≈ code_phase

    @test signal ≈ cis.(2π * (interm_freq + doppler) / sample_freq * test_range .+ carrier_phase) .* gen_code.(Ref(gnss_system), test_range, code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0) * 1/1Hz)

    @test @inferred(GNSSSimulator.calc_code_phase(0m, 1_023_000Hz, 1023)) == 0
    @test @inferred(GNSSSimulator.calc_code_phase(293.255132m, 1_023_000Hz, 1023)) ≈ 1 rtol = 1e-3 # 1 Chip is around 300m

    @test @inferred(GNSSSimulator.calc_amplitude_from_cn0(45dBHz, 1/1Hz)) ≈ 10^(45 / 20)

    next_measurement, signal, internal_states = @inferred next_measurement(num_samples)
    @test signal ≈ cis.(2π .* (interm_freq .+ doppler) ./ sample_freq .* (num_samples + 1:2 * num_samples) .+ carrier_phase) .* gen_code.(Ref(gnss_system), (num_samples + 1:2 * num_samples), code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0) * 1/1Hz)
end

@testset "Satellite L5" begin
    cn0 = 10.0dBHz
    sample_freq = 40e6Hz
    interm_freq = 100Hz
    center_freq = 1_176_450_000Hz
    code_length = 102300
    num_samples = 400000
    test_range = 1:num_samples
    gnss_system = GPSL5()
    attitude = RotXYZ(0.0, 0.0, 0.0)
    sat = Satellite(prn = 1, enu_doa = CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), CN0 = cn0)
    measurement, init_internal_states = @inferred GNSSSimulator.init_sim_emitter_signal(sat, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(num_samples)
    doppler = calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, center_freq)
    code_freq = 10.23MHz + doppler * 10.23MHz / center_freq
    sat_user_distance = calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    code_phase = GNSSSimulator.calc_code_phase(sat_user_distance, code_freq, code_length)
    carrier_phase = GNSSSimulator.calc_carrier_phase(sat_user_distance, center_freq + doppler)

    @test init_internal_states.doppler ≈ doppler
    @test init_internal_states.carrier_phase ≈ carrier_phase
    @test init_internal_states.code_phase ≈ code_phase

    @test signal ≈ cis.(2π * (interm_freq + doppler) / sample_freq * test_range .+ carrier_phase) .* gen_code.(Ref(gnss_system), test_range, code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0) * 1/1Hz)
end
