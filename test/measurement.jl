@testset "Measurement" begin

    @testset "Emitters" begin
        manifold = IdealManifold()
        sample_freq = 2e6Hz

        sat1 = ConstantDopplerSatellite(GPSL1, 1, carrier_doppler = 1000.0Hz, carrier_phase = π / 2, code_phase = 100.0, cn0 = 45dBHz)
        sat2 = ConstantDopplerSatellite(GPSL1, 2, carrier_doppler = 500.0Hz, carrier_phase = π / 2, code_phase = 50.0, cn0 = 45dBHz)
        emitters = (sat1,sat2)
        receiver = @inferred Receiver(1.0, RotXYZ(0.0, 0.0, 0.0), 0.0)
        received_signal = @inferred ReceivedSignal(emitters, receiver)

        next_received_signal = @inferred propagate(received_signal, 1 / sample_freq)

        next_sat = next_received_signal.emitters[1]
        @test GNSSSimulator.get_carrier_doppler(next_sat) == 1000Hz
        @test GNSSSimulator.get_code_doppler(next_sat) ≈ 1000Hz / 1540
        @test GNSSSimulator.get_carrier_phase(next_sat) ≈ π / 2 + 2π * 1000Hz / 2e6Hz
        @test GNSSSimulator.get_code_phase(next_sat) ≈ 100 + (1023e3Hz + 1000Hz / 1540) / 2e6Hz
        @test GNSSSimulator.get_existence(next_sat) == true
        @test GNSSSimulator.get_amplitude(next_sat) == 10^(45 / 20)
        @test GNSSSimulator.get_prn(next_sat) == 1

        measurement = @inferred get_measurement(next_received_signal, manifold)
        @test measurement ≈ cis(π / 2 + 2π * 1000Hz / 2e6Hz) * get_code(GPSL1, 100 + (1023e3Hz + 1000Hz / 1540) / 2e6Hz, 1) * 10^(45 / 20) +
            cis(π / 2 + 2π * 500Hz / 2e6Hz) * get_code(GPSL1, 50 + (1023e3Hz + 500Hz / 1540) / 2e6Hz, 2) * 10^(45 / 20)
    end

    @testset "Noise" begin
        manifold = IdealManifold()
        bandwidth = 2e6Hz
        sat = ConstantDopplerSatellite(GPSL1, 1, cn0 = 0.0dBHz) # Dummy satellite
        emitters = (sat,)
        receiver = Receiver(1.0, RotXYZ(0.0, 0.0, 0.0), sqrt(bandwidth * 1.0/Hz))
        received_signal = ReceivedSignal(emitters, receiver)
        measurements = [get_measurement(propagate(received_signal, 1μs), manifold) for i = 1:1000]
        @test measurements'measurements / 1000 ≈ 2e6 atol = 0.03e6
    end
end
