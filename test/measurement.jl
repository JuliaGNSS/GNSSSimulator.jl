@testset "Continuity of measurement" begin 
    sample_freq = 4e6
    interm_freq = 100_000
    gnss_system = GNSSSimulator.gpsl1_system()
    jammer = GNSSSimulator.CWJammer(1, Spherical(0.0, 0.0, 1.0), 0.0, 20.0)
    sat = GNSSSimulator.Satellite(svid = 1, enu_doa = Spherical(0.0, 0.0, 1.0))
    attitude = RotXYZ(0.0, 0.0, 0.0)
    emitters = [sat, jammer]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal1, internal_states = measurement(10)
    next_measurement, signal_short1, internal_states = measurement(5)
    next_measurement, signal_short2, internal_states = next_measurement(5)
    signal2 = [signal_short1 signal_short2]
    @test signal1 ≈ signal2
end

@testset "Measurement" begin 
    sample_freq = 4e6
    interm_freq = 100_000
    jammer_power_dB = 20.0
    gnss_system = GNSSSimulator.gpsl1_system()
    jammer = GNSSSimulator.CWJammer(1, Spherical(0.0, 0.0, 1.0), 0.0, jammer_power_dB)
    attitude = RotXYZ(0.0, 0.0, 0.0)
    emitters = [jammer]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement = @inferred GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq, false)

    next_measurement, signal, internal_states = measurement(100)
    @test signal ≈ get_steer_vec(Spherical(0.0, 0.0, 1.0)) * cis.(2π * interm_freq / sample_freq * (1:100)') * 10^(jammer_power_dB / 20)
end