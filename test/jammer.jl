@testset "Jammer" begin
    jammer_power_dB = 10.0
    sample_freq = 4e6
    interm_freq = 100_000
    gnss_system = GNSSSimulator.gpsl1_system()
    attitude = RotXYZ(0.0, 0.0, 0.0)
    jammer = CWJammer(1, Spherical(0.0, 0.0, 1.0), 0.0, jammer_power_dB)
    measurement = @inferred GNSSSimulator.init_sim_emitter_signal(jammer, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(100)
    @test signal ≈ cis.(2π * interm_freq / sample_freq * (1:100)) * 10^(jammer_power_dB / 20)

    @test GNSSSimulator.calc_amplitude_from_jnr(10) ≈ 10^(10 / 20)
end