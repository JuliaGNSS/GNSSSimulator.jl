@testset "Structural interference" begin
    gpsl1 = GPSL1()
    sat = ConstantDopplerSatellite(GPSL1, 1, carrier_doppler = 1000.0Hz, carrier_phase = π / 2, code_phase = 100.0, cn0 = 45dBHz)
    si = ConstantDopplerStructuralInterference(sat, -3dB, added_carrier_doppler = 10Hz, added_carrier_phase = π / 8, added_code_phase = 10.0)

    phase = @inferred GNSSSimulator.get_phase(si)
    @test phase.carrier == π / 2 + π / 8
    @test phase.code == 110.0

    signal = @inferred GNSSSimulator.get_signal(phase, si, 1.0 + 0.0im, Random.GLOBAL_RNG)
    @test signal ≈ 1.0 * 10^(42 / 20) * GNSSSignals.cis_vfast(π / 2 + π / 8) * get_code(GPSL1, 110.0, 1)

    signal = @inferred GNSSSimulator.get_signal(phase, si, SVector(1.0im, 2.0im), Random.GLOBAL_RNG)
    @test signal ≈ SVector(1im, 2im) * 10^(42 / 20) * GNSSSignals.cis_vfast(π / 2 + π / 8) * get_code(GPSL1, 110.0, 1)

    next_phase = @inferred GNSSSimulator.fast_propagate(phase, si, 100.0Hz, 1/2e6Hz)
    @test next_phase.carrier == π / 2 + π / 8 + 2π * 1110.0Hz / 2e6Hz
    @test next_phase.code == 110.0 + (1023e3Hz + 1010.0Hz / 1540) / 2e6Hz

    next_si = @inferred GNSSSimulator.propagate(si, 100.0Hz, 1/2e6Hz, Random.GLOBAL_RNG)
    @test @inferred(get_carrier_doppler(next_si)) == 1010Hz
    @test @inferred(get_code_doppler(next_si)) ≈ 1010Hz / 1540
    @test @inferred(get_carrier_phase(next_si)) ≈ π / 2 + π / 8 + 2π * 1110Hz / 2e6Hz
    @test @inferred(get_code_phase(next_si)) ≈ 110 + (1023e3Hz + 1010Hz / 1540) / 2e6Hz
    @test @inferred(get_existence(next_si)) == true
    @test @inferred(get_amplitude(next_si)) ≈ 10^(42 / 20)
    @test @inferred(get_prn(next_si)) == 1
end
