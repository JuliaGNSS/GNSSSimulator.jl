@testset "Jammer" begin

    @testset "Phase wrap" begin
        jammer = CWJammer(1, 10)
        phase = GNSSSimulator.CWJammerPhase(0.25)
        phase_wrap = GNSSSimulator.CWJammerPhaseWrap(0)

        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(
            phase_wrap,
            phase,
            jammer
        )
        @test next_phase_wrap == GNSSSimulator.CWJammerPhaseWrap(0)

        phase = GNSSSimulator.CWJammerPhase(0.55)
        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(
            next_phase_wrap,
            phase,
            jammer
        )
        @test next_phase_wrap == GNSSSimulator.CWJammerPhaseWrap(1)
    end

    @testset "CW Jammer" begin
        jammer = CWJammer(1, 100.0Hz, 0.5, 10.0, true, SVector(0.0, 0.0, 1.0))

        @test @inferred(GNSSSimulator.get_carrier_phase(jammer)) ≈ π

        phase = @inferred GNSSSimulator.calc_phase(jammer, 0, 100.0Hz, 1e6Hz)
        @test phase.carrier ≈ 0.5

        phase_wrap = @inferred GNSSSimulator.init_phase_wrap(jammer)
        @test phase_wrap == GNSSSimulator.CWJammerPhaseWrap(0)

        signal = @inferred GNSSSimulator.get_signal(
            jammer,
            phase,
            phase_wrap,
            1.0 + 0.0im,
            Random.GLOBAL_RNG
        )
        @test signal ≈ 1.0 * 10 * GNSSSignals.cis_vfast(1.0π)

        next_phase = @inferred GNSSSimulator.calc_phase(jammer, 1, 100.0Hz, 1e6Hz)
        @test next_phase.carrier ≈ 0.5 + 200.0Hz * 1μs

        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(
            phase_wrap,
            next_phase,
            jammer
        )

        @test next_phase_wrap == GNSSSimulator.CWJammerPhaseWrap(1)

        signal = @inferred GNSSSimulator.get_signal(
            jammer,
            next_phase,
            next_phase_wrap,
            1.0 + 0.0im,
            Random.GLOBAL_RNG
        )
        @test signal ≈ 1.0 * 10 * GNSSSignals.cis_vfast(π + 2π * 200.0Hz * 1μs) rtol = 1

        next_jammer = @inferred GNSSSimulator.propagate(
            jammer,
            100.0Hz,
            1μs,
            Random.GLOBAL_RNG
        )
        @test @inferred(get_carrier_phase(next_jammer)) ≈ mod2pi(
            π + 2π * 200.0Hz * 1μs + π
        ) - π
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

        phase = @inferred GNSSSimulator.calc_phase(jammer, 0, 100.0Hz, 1e6Hz)
        @test phase == GNSSSimulator.NoiseJammerPhase()

        phase_wrap = @inferred GNSSSimulator.init_phase_wrap(jammer)
        @test phase_wrap == GNSSSimulator.NoiseJammerPhaseWrap()

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.get_signal(
            jammer,
            phase,
            phase_wrap,
            SVector(1.0im, 2.0im),
            rng
        )
        rng = MersenneTwister(1234)
        @test signal ≈ SVector(1im, 2im) * 10 * randn(rng, ComplexF64)

        next_phase = @inferred GNSSSimulator.calc_phase(jammer, 1, 100.0Hz, 1e6Hz)
        @test next_phase == GNSSSimulator.NoiseJammerPhase()

        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(
            phase_wrap,
            next_phase,
            jammer
        )
        @test next_phase_wrap == GNSSSimulator.NoiseJammerPhaseWrap()

        next_jammer = @inferred GNSSSimulator.propagate(
            jammer,
            100.0Hz,
            1μs,
            Random.GLOBAL_RNG
        )
        @test @inferred(get_amplitude(next_jammer)) == 10
        @test @inferred(get_existence(next_jammer)) == true
        @test @inferred(get_doa(next_jammer)) == SVector(0.0, 0.0, 1.0)
        @test @inferred(get_id(next_jammer)) == 1

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.get_signal(
            jammer,
            next_phase,
            next_phase_wrap,
            SVector(1.0im, 2.0im),
            rng
        )
        rng = MersenneTwister(1234)
        @test signal ≈ SVector(1im, 2im) * 10 * randn(rng, ComplexF64)
    end
end
