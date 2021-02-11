@testset "Measurement" begin

    @testset "Multiple Emitters" begin
        manifold = IdealManifold()
        system = GPSL1()

        sat1 = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sat2 = ConstantDopplerSatellite(
            system,
            2,
            carrier_doppler = 500.0Hz,
            carrier_phase = π / 2,
            code_phase = 50.0,
            cn0 = 45dBHz
        )
        sats = (sat1,sat2)
        receiver = @inferred Receiver(2.5e6Hz)

        signal = StructArray{ComplexF32}(undef, 2500)
        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement!(
            signal,
            receiver,
            sats,
            manifold,
            rng
        )

        carrier = cis.(2π .* 1000.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal = carrier .*
            get_code.(system, (0:2499) .* (1023e3 + 1000 / 1540) ./ 2.5e6 .+ 100, 1) .*
            sqrt(10^(45 / 10))
        carrier = cis.(2π .* 500.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal += carrier .*
            get_code.(system, (0:2499) .* (1023e3 + 500 / 1540) ./ 2.5e6 .+ 50, 2) .*
            sqrt(10^(45 / 10))
        rng = MersenneTwister(1234)
        reference_signal += randn(rng, ComplexF64, 2500) .* sqrt(2.5e6)
        @test measurement ≈ reference_signal

        @test next_receiver == GNSSSimulator.propagate(receiver, 2500, rng)
        @test next_sats == GNSSSimulator.propagate.(sats, 2500, 0.0Hz, 2.5e6Hz, Ref(rng))

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement(
            2500,
            receiver,
            sats,
            manifold,
            rng
        )
        @test measurement ≈ reference_signal

        sat1f64 = ConstantDopplerSatellite(
            Float64,
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sat2f64 = ConstantDopplerSatellite(
            Float64,
            system,
            2,
            carrier_doppler = 500.0Hz,
            carrier_phase = π / 2,
            code_phase = 50.0,
            cn0 = 45dBHz
        )
        satsf64 = (sat1f64,sat2f64)
        receiver = @inferred Receiver(2.5e6Hz)

        signalf64 = StructArray{ComplexF64}(undef, 2500)
        rng = MersenneTwister(1234)
        measurementf64, next_receiverf64, next_satsf64 = @inferred get_measurement!(
            signalf64,
            receiver,
            satsf64,
            manifold,
            rng
        )
        @test measurementf64 ≈ reference_signal

        rng = MersenneTwister(1234)
        measurementf64, next_receiver, next_sats = @inferred get_measurement(
            2500,
            receiver,
            satsf64,
            manifold,
            rng
        )
        @test measurementf64 ≈ reference_signal
    end

    @testset "Non existing Emitter" begin
        manifold = IdealManifold()
        system = GPSL1()

        sat1 = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sat2 = ConstantDopplerSatellite(
            system,
            2,
            carrier_doppler = 500.0Hz,
            carrier_phase = π / 2,
            code_phase = 50.0,
            cn0 = 45dBHz,
            exists = false
        )
        sats = (sat1,sat2)
        receiver = @inferred Receiver(2.5e6Hz)

        signal = StructArray{ComplexF32}(undef, 2500)
        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement!(
            signal,
            receiver,
            sats,
            manifold,
            rng
        )

        carrier = cis.(2π .* 1000.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal = carrier .*
            get_code.(system, (0:2499) .* (1023e3 + 1000 / 1540) ./ 2.5e6 .+ 100, 1) .*
            sqrt(10^(45 / 10))
        rng = MersenneTwister(1234)
        reference_signal += randn(rng, ComplexF32, 2500) .* sqrt(2.5e6)
        @test measurement ≈ reference_signal
    end

    @testset "Signal continuity" begin
        manifold = IdealManifold()
        system = GPSL1()

        sat1 = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1,)
        receiver = @inferred Receiver(2.5e6Hz)

        rng = MersenneTwister(1234)
        measurement1, next_receiver, next_sats = @inferred get_measurement(
            2500,
            receiver,
            sats,
            manifold,
            rng
        )

        @test next_receiver == GNSSSimulator.propagate(receiver, 2500, rng)
        @test next_sats == GNSSSimulator.propagate.(sats, 2500, 0.0Hz, 2.5e6Hz, Ref(rng))

        measurement2, next_receiver, next_sats = @inferred get_measurement(
            2500,
            next_receiver,
            next_sats,
            manifold,
            rng
        )
        measurement = [measurement1; measurement2]

        carrier = cis.(2π .* 1000.0Hz ./ 2.5e6Hz .* (0:4999) .+ 2π .* 0.25)
        reference_signal = carrier .*
            get_code.(system, (0:4999) .* (1023e3 + 1000 / 1540) ./ 2.5e6 .+ 100, 1) .*
            sqrt(10^(45 / 10))
        rng = MersenneTwister(1234)
        reference_signal += randn(rng, ComplexF32, 5000) .* sqrt(2.5e6)

        @test measurement.im ≈ imag.(reference_signal)
        # Slightly not contineous because of Float and Fixed point representation?
        @test measurement.re[1:2500] ≈ real.(reference_signal[1:2500])
        @test measurement.re[2502:4917] ≈ real.(reference_signal[2502:4917])
        @test measurement.re[4919:end] ≈ real.(reference_signal[4919:end])
    end

    @testset "Multiple antennas" begin
        manifold = IdealManifold(
            0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0)),
            1575420e3
        )
        system = GPSL1()

        sat1 = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1,)
        gpmc = @SMatrix randn(ComplexF64, 2, 2)
        receiver = @inferred Receiver(
            2.5e6Hz,
            gain_phase_mism_crosstalk = gpmc
        )

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement(
            2500,
            receiver,
            sats,
            manifold,
            rng
        )
        carrier = cis.(2π .* 1000.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal = carrier .*
            get_code.(system, (0:2499) .* (1023e3 + 1000 / 1540) ./ 2.5e6 .+ 100, 1) .*
            sqrt(10^(45 / 10)) .* transpose([1.0, 1.0])
        rng = MersenneTwister(1234)
        reference_signal += randn(rng, ComplexF32, 2500, 2) .* sqrt(2.5e6)
        reference_signal = reference_signal * transpose(gpmc)
        @test measurement ≈ reference_signal
    end

    @testset "Two emitter types" begin
        manifold = IdealManifold()
        gpsl1 = GPSL1()
        gpsl5 = GPSL5()

        sat1 = ConstantDopplerSatellite(
            gpsl1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sat2 = ConstantDopplerSatellite(
            gpsl5,
            1,
            carrier_doppler = 500.0Hz,
            carrier_phase = π / 2,
            code_phase = 50.0,
            cn0 = 45dBHz
        )
        sats = (sat1,sat2)
        receiver = @inferred Receiver(2.5e6Hz)

        jammer = CWJammer(1, 10dB, doppler = 100.0Hz, phase = π / 4)
        jammers = (jammer, )

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats, next_jammers = @inferred get_measurement(
            2500,
            receiver,
            sats,
            jammers,
            manifold,
            rng
        )

        carrier = cis.(2π .* 1000.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal = carrier .*
            get_code.(gpsl1, (0:2499) .* (1023e3 + 1000 / 1540) ./ 2.5e6 .+ 100, 1) .*
            sqrt(10^(45 / 10))
        carrier = cis.(2π .* 500.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal += carrier .*
            get_code.(gpsl5, (0:2499) .* (1023e4 + 500 / 115) ./ 2.5e6 .+ 50, 1) .*
            sqrt(10^(45 / 10))
        rng = MersenneTwister(1234)
        reference_signal += randn(rng, ComplexF64, 2500) .* sqrt(2.5e6)
        carrier = cis.(2π .* 100.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.125)
        reference_signal += carrier .* sqrt(10^(10 / 10) * 2.5e6)
        @test measurement ≈ reference_signal

    end

    @testset "Three emitter types" begin
        manifold = IdealManifold()
        system = GPSL1()

        sat1 = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1, )
        receiver = @inferred Receiver(2.5e6Hz)

        jammer = CWJammer(1, 10dB, doppler = 100.0Hz, phase = π / 4)
        jammers = (jammer, )

        si = ConstantDopplerStructuralInterference(
            sat1,
            -3dB,
            added_carrier_doppler = 10Hz,
            added_carrier_phase = π / 8,
            added_code_phase = 10.0
        )

        sis = (si, )

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats, next_jammers, next_sis = @inferred get_measurement(
            2500,
            receiver,
            sats,
            jammers,
            sis,
            manifold,
            rng
        )

        carrier = cis.(2π .* 1000.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.25)
        reference_signal = carrier .*
            get_code.(system, (0:2499) .* (1023e3 + 1000 / 1540) ./ 2.5e6 .+ 100, 1) .*
            sqrt(10^(45 / 10))
        carrier = cis.(2π .* 1010.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.3125)
        reference_signal += carrier .*
            get_code.(system, (0:2499) .* (1023e3 + 1010 / 1540) ./ 2.5e6 .+ 110, 1) .*
            sqrt(10^(42 / 10))
        rng = MersenneTwister(1234)
        reference_signal += randn(rng, ComplexF64, 2500) .* sqrt(2.5e6)
        carrier = cis.(2π .* 100.0Hz ./ 2.5e6Hz .* (0:2499) .+ 2π .* 0.125)
        reference_signal += carrier .* sqrt(10^(10 / 10) * 2.5e6)
        @test measurement ≈ reference_signal
    end
end
