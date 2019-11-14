@testset "Measurement" begin

    @testset "Multiple Emitters" begin
        manifold = IdealManifold()

        sat1 = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sat2 = ConstantDopplerSatellite(
            GPSL1,
            2,
            carrier_doppler = 500.0Hz,
            carrier_phase = π / 2,
            code_phase = 50.0,
            cn0 = 45dBHz
        )
        sats = (sat1,sat2)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)

        signal = Vector{ComplexF64}(undef, 4000)
        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement!(
            signal,
            receiver,
            sats,
            manifold,
            rng
        )

        x = 0:3999
        rng = MersenneTwister(1234)
        test_signal =
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ 2π * 1000Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 100 .+ (1023e3Hz + 1000Hz / 1540) / 4e6Hz * x, 1) .*
            10^(45 / 20) .+
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ 2π * 500Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 50 .+ (1023e3Hz + 500Hz / 1540) / 4e6Hz * x, 2) .*
            10^(45 / 20) .+
            randn(rng, ComplexF64, 4000) .* 5.0

        @test measurement ≈ test_signal

        @test next_receiver == GNSSSimulator.propagate(receiver, 4000, rng)
        @test next_sats == GNSSSimulator.propagate(sats, 4000, 0.0Hz, 4e6Hz, rng)

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement(
            4000,
            receiver,
            sats,
            manifold,
            rng
        )
        @test measurement ≈ test_signal

        rng = MersenneTwister(1234)
        measurement_f32, next_receiver, next_sats = @inferred get_measurement(
            Float32,
            4000,
            receiver,
            sats,
            manifold,
            rng
        )
        @test measurement_f32 ≈ test_signal

        rng = MersenneTwister(1234)
        @test next_receiver == GNSSSimulator.propagate(receiver, 4000, rng)
        @test next_sats == GNSSSimulator.propagate(sats, 4000, 0.0Hz, 4e6Hz, rng)
    end

    @testset "Non existing Emitter" begin
        manifold = IdealManifold()
        sat1 = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sat2 = ConstantDopplerSatellite(
            GPSL1,
            2,
            carrier_doppler = 500.0Hz,
            carrier_phase = π / 2,
            code_phase = 50.0,
            cn0 = 45dBHz,
            exists = false
        )
        sats = (sat1,sat2)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement(
            4000,
            receiver,
            sats,
            manifold,
            rng
        )

        x = 0:3999
        rng = MersenneTwister(1234)
        test_signal =
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ 2π * 1000Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 100 .+ (1023e3Hz + 1000Hz / 1540) / 4e6Hz * x, 1) .*
            10^(45 / 20) .+
            randn(rng, ComplexF64, 4000) .* 5.0
        @test measurement ≈ test_signal
    end

    @testset "Signal continuity" begin
        manifold = IdealManifold()
        sat1 = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1,)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)

        rng = MersenneTwister(1234)
        measurement1, next_receiver, next_sats = @inferred get_measurement(
            4000,
            receiver,
            sats,
            manifold,
            rng
        )

        @test next_receiver == GNSSSimulator.propagate(receiver, 4000, rng)
        @test next_sats == GNSSSimulator.propagate(sats, 4000, 0.0Hz, 4e6Hz, rng)

        measurement2, next_receiver, next_sats = @inferred get_measurement(
            4000,
            next_receiver,
            next_sats,
            manifold,
            rng
        )
        measurement = [measurement1; measurement2]

        x = 0:7999
        rng = MersenneTwister(1234)
        test_signal =
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ 2π * 1000Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 100 .+ (1023e3Hz + 1000Hz / 1540) / 4e6Hz * x, 1) .*
            10^(45 / 20) .+
            randn(rng, ComplexF64, 8000) .* 5.0
        @test measurement ≈ test_signal
    end

    @testset "Multiple antennas" begin
        manifold = IdealManifold(
            0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0)),
            1575420e3
        )
        sat1 = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1,)
        receiver = @inferred Receiver(
            4e6Hz,
            noise_std = 5.0,
            gain_phase_mism_crosstalk = SMatrix{2,2,ComplexF64}(I)
        )

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats = @inferred get_measurement(
            4000,
            receiver,
            sats,
            manifold,
            rng
        )
        x = 0:3999
        rng = MersenneTwister(1234)
        test_signal =
            [1.0; 1.0] .* transpose(
                GNSSSignals.cis_vfast.(
                    mod2pi.(π / 2 .+ 2π * 1000Hz / 4e6Hz * x .+ π) .- π
                ) .*
                get_code.(GPSL1, 100 .+ (1023e3Hz + 1000Hz / 1540) / 4e6Hz * x, 1) .*
                10^(45 / 20)
            ) .+ randn(rng, ComplexF64, 2, 4000) .* 5.0
        @test measurement ≈ test_signal
    end

    @testset "Two emitter types" begin

        manifold = IdealManifold()
        sat1 = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1,)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)

        jammer = CWJammer(1, 10, doppler = 100.0Hz)
        jammers = (jammer, )

        rng = MersenneTwister(1234)
        measurement, next_receiver, next_sats, next_jammers = @inferred get_measurement(
            4000,
            receiver,
            sats,
            jammers,
            manifold,
            rng
        )

        x = 0:3999
        rng = MersenneTwister(1234)
        test_signal =
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ 2π * 1000Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 100 .+ (1023e3Hz + 1000Hz / 1540) / 4e6Hz * x, 1) .*
            10^(45 / 20) .+ 10 .*
            GNSSSignals.cis_vfast.(mod2pi.(2π * 100.0Hz / 4e6Hz * x .+ π) .- π) .+
            randn(rng, ComplexF64, 4000) .* 5.0

        @test measurement ≈ test_signal
    end

    @testset "Three emitter types" begin

        manifold = IdealManifold()
        sat1 = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )
        sats = (sat1,)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)

        jammer = CWJammer(1, 10, doppler = 100.0Hz)
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
            4000,
            receiver,
            sats,
            sis,
            jammers,
            manifold,
            rng
        )

        x = 0:3999
        rng = MersenneTwister(1234)
        test_signal =
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ 2π * 1000Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 100 .+ (1023e3Hz + 1000Hz / 1540) / 4e6Hz * x, 1) .*
            10^(45 / 20) .+
            GNSSSignals.cis_vfast.(mod2pi.(π / 2 .+ π / 8 .+ 2π * 1010Hz / 4e6Hz * x .+ π) .- π) .*
            get_code.(GPSL1, 110 .+ (1023e3Hz + 1010Hz / 1540) / 4e6Hz * x, 1) .*
            10^(42 / 20) .+ 10 .*
            GNSSSignals.cis_vfast.(mod2pi.(2π * 100.0Hz / 4e6Hz * x .+ π) .- π) .+
            randn(rng, ComplexF64, 4000) .* 5.0

        @test measurement ≈ test_signal
    end


end
