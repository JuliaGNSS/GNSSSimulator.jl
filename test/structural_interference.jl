@testset "Structural interference" begin
    gpsl1 = GPSL1()
    sat = ConstantDopplerSatellite(1, gpsl1, carrier_doppler = 1000.0Hz, carrier_phase = π / 2, code_phase = 100.0, cn0 = 45dBHz)
    si = ConstantDopplerStructuralInterference(sat, -3dB, added_carrier_doppler = 10Hz, added_carrier_phase = π / 8, added_code_phase = 10)
    next_si = @inferred propagate(si, 1/2e6Hz)
    @test @inferred(GNSSSimulator.get_carrier_doppler(next_si)) == 1010Hz
    @test @inferred(GNSSSimulator.get_code_doppler(next_si)) ≈ 1010Hz / 1540
    @test @inferred(GNSSSimulator.get_carrier_phase(next_si)) ≈ π / 2 + π / 8 + 2π * 1010Hz / 2e6Hz
    @test @inferred(GNSSSimulator.get_code_phase(next_si)) ≈ 110 + 1010Hz / 1540 / 2e6Hz
    @test @inferred(GNSSSimulator.get_existence(next_si)) == true
    @test @inferred(GNSSSimulator.get_amplitude(next_si)) ≈ 10^(42 / 20)
    @test @inferred(GNSSSimulator.get_prn(next_si)) == 1
    @test @inferred(GNSSSimulator.get_system(next_si)) == gpsl1

    get_steer_vec(doa, attitude) = 1 # No steering vector

    signal = @inferred get_signal(next_si, nothing, get_steer_vec)
    @test signal ≈ cis(π / 2 + π / 8 + 2π * 1010Hz / 2e6Hz) * gpsl1.codes[1 + floor(Int, 110 + 1010Hz / 1540 / 2e6Hz), 1] * 10^(42 / 20)
end
