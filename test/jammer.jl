@testset "Jammer" begin
    @testset "CW Jammer" begin
        jammer = CWJammer(1, 100.0Hz, 1.0, 10.0, true, SVector(0.0, 0.0, 1.0))

        phase = @inferred GNSSSimulator.get_phase(jammer)
        @test phase.carrier == 1.0

        signal = @inferred GNSSSimulator.get_signal(phase, jammer, 1.0 + 0.0im, Random.GLOBAL_RNG)
        @test signal ≈ 1.0 * 10 * GNSSSignals.cis_vfast(1.0)

        next_phase = @inferred GNSSSimulator.fast_propagate(phase, jammer, 100.0Hz, 1μs)
        @test next_phase.carrier == 1.0 + 2π * 200.0Hz * 1μs

        signal = @inferred GNSSSimulator.get_signal(next_phase, jammer, 1.0 + 0.0im, Random.GLOBAL_RNG)
        @test signal ≈ 1.0 * 10 * GNSSSignals.cis_vfast(1.0 + 2π * 200.0Hz * 1μs)

        next_jammer = @inferred GNSSSimulator.propagate(jammer, 100.0Hz, 1μs, Random.GLOBAL_RNG)
        @test @inferred(get_carrier_phase(next_jammer)) ≈ 1.0 + 2π * 200.0Hz * 1μs
        @test @inferred(get_carrier_doppler(next_jammer)) == 100.0Hz
        @test @inferred(get_amplitude(next_jammer)) == 10
        @test @inferred(get_existence(next_jammer)) == true
        @test @inferred(get_id(next_jammer)) == 1
        @test @inferred(get_doa(next_jammer)) == SVector(0.0, 0.0, 1.0)

        jammer = @inferred CWJammer(1, 10.0)
        @test get_id(jammer) == 1
        @test get_amplitude(jammer) == 10.0
    end

    @testset "Noise Jammer" begin
        jammer = @inferred NoiseJammer(1, 10.0, true, SVector(0.0, 0.0, 1.0))

        phase = @inferred GNSSSimulator.get_phase(jammer)
        @test phase == GNSSSimulator.NoiseJammerPhase()

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.get_signal(phase, jammer, SVector(1.0im, 2.0im), rng)
        rng = MersenneTwister(1234)
        @test signal ≈ SVector(1im, 2im) * 10 * randn(rng, ComplexF64)

        next_phase = @inferred GNSSSimulator.fast_propagate(phase, jammer, 100.0Hz, 1μs)
        @test next_phase == GNSSSimulator.NoiseJammerPhase()

        next_jammer = @inferred GNSSSimulator.propagate(jammer, 100.0Hz, 1μs, Random.GLOBAL_RNG)
        @test @inferred(get_amplitude(next_jammer)) == 10
        @test @inferred(get_existence(next_jammer)) == true
        @test @inferred(get_doa(next_jammer)) == SVector(0.0, 0.0, 1.0)
        @test @inferred(get_id(next_jammer)) == 1

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.get_signal(next_phase, jammer, SVector(1.0im, 2.0im), rng)
        rng = MersenneTwister(1234)
        @test signal ≈ SVector(1im, 2im) * 10 * randn(rng, ComplexF64)
    end
end
