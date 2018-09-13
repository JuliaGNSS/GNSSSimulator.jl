@testset "Jammer" begin
    jammer_power = 10.0dB
    sample_freq = 4e6Hz
    interm_freq = 100_000Hz
    gnss_system = GPSL1()
    jammer = CWJammer(1, CartesianFromSpherical()(Spherical(0.0, 0.0, 1.0)), 0.0m / 1s, jammer_power, true)
    measurement = @inferred GNSSSimulator.init_sim_emitter_signal(jammer, gnss_system, sample_freq, interm_freq)
    next_measurement, signal, internal_states = @inferred measurement(100)
    @test signal ≈ cis.(2π * interm_freq / sample_freq * (1:100)) * sqrt(uconvertp(NoUnits, jammer_power) * sample_freq / 1Hz)
end