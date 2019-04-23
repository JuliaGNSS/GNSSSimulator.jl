@testset "Structural interference for satellite L1" begin
    cn0 = 45.0dBHz
    sample_freq = 4e6Hz
    interm_freq = 100_000Hz
    center_freq = 1_575_420_000Hz
    code_length = 1023
    num_samples = 4000
    added_relative_velocity = 5.0m/s
    signal_amplification = -3dB
    added_signal_path = 10m
    test_range = 1:num_samples
    gnss_system = GPSL1()
    attitude = RotXYZ(0.0, 0.0, 0.0)
    sat = Satellite(prn = 1, enu_doa = CartesianFromSpherical()(Spherical(1.0, 0.0, 0.0)), CN0 = cn0)
    structural_interference = StructuralInterference(enu_doa = CartesianFromSpherical()(Spherical(1.0, π/2, 0.0)), sat = sat, signal_amplification = signal_amplification, added_relative_velocity = added_relative_velocity, added_signal_path = added_signal_path)
    measurement, init_internal_states = @inferred GNSSSimulator.init_sim_emitter_signal(structural_interference, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(num_samples)
    doppler = calc_init_doppler(sat.distance_from_earth_center, sat.enu_doa, sat.velocity, center_freq) + added_relative_velocity / GNSSSimulator.SPEED_OF_LIGHT * gnss_system.center_freq
    code_freq = 1_023_000Hz + doppler * 1_023_000Hz / center_freq
    sat_user_distance = calc_init_sat_user_distance(sat.distance_from_earth_center, sat.enu_doa)
    code_phase = GNSSSimulator.calc_code_phase(sat_user_distance + added_signal_path, code_freq, code_length)
    carrier_phase = GNSSSimulator.calc_carrier_phase(sat_user_distance + added_signal_path, center_freq + doppler)

    @test signal' * signal / 4000 ≈ 10^(42 / 10)

    @test init_internal_states.doppler ≈ doppler
    @test init_internal_states.carrier_phase ≈ carrier_phase
    @test init_internal_states.code_phase ≈ code_phase

    @test signal ≈ cis.(2π * (interm_freq + doppler) / sample_freq * test_range .+ carrier_phase) .* gen_code.(Ref(gnss_system), test_range, code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0 + signal_amplification) * 1/1Hz)

    next_measurement, signal, internal_states = @inferred next_measurement(num_samples)
    @test signal ≈ cis.(2π .* (interm_freq .+ doppler) ./ sample_freq .* (num_samples + 1:2 * num_samples) .+ carrier_phase) .* gen_code.(Ref(gnss_system), (num_samples + 1:2 * num_samples), code_freq, code_phase, sample_freq, 1) .* sqrt(linear(cn0 + signal_amplification) * 1/1Hz)
end
