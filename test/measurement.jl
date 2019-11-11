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
        emitters = (sat1,sat2)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)

        signal = Vector{ComplexF64}(undef, 4000)
        rng = MersenneTwister(1234)
        measurement = @inferred get_measurement!(signal, receiver, emitters, manifold, rng)

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

        using PyPlot
        pygui(true)
        plot(real.(measurement))
        plot(real.(test_signal))

        @test measurement ≈ test_signal

        rng = MersenneTwister(1234)
        @test get_measurement(4000, receiver, emitters, manifold, rng) ≈ test_signal
        rng = MersenneTwister(1234)
        @test get_measurement(
            Float32,
            4000,
            receiver,
            emitters,
            manifold,
            rng
        ) ≈ test_signal

        next_receiver = propagate(receiver, 1ms)
        next_emitters = propagate(emitters, 0.0Hz, 1ms)
        next_sat1 = next_emitters[1]

        @test get_carrier_doppler(next_sat1) == 1000Hz
        @test get_code_doppler(next_sat1) ≈ 1000Hz / 1540
        @test get_carrier_phase(next_sat1) ≈ mod2pi(
            π / 2 + 2π * 1000Hz / 4e6Hz * 4000 + π
        ) - π
        @test get_code_phase(next_sat1) ≈ mod(
            100 + (1023e3Hz + 1000Hz / 1540) / 4e6Hz * 4000,
            1023
        )
        @test get_existence(next_sat1) == true
        @test get_amplitude(next_sat1) == 10^(45 / 20)
        @test get_prn(next_sat1) == 1
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
        emitters = (sat1,sat2)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)
        received_signal = @inferred ReceivedSignal(emitters, receiver)

        rng = MersenneTwister(1234)
        measurement = @inferred get_measurement(4000, received_signal, manifold, rng)

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
        emitters = (sat1,)
        receiver = @inferred Receiver(4e6Hz, noise_std = 5.0)
        received_signal = @inferred ReceivedSignal(emitters, receiver)

        rng = MersenneTwister(1234)
        measurement1 = @inferred get_measurement(4000, received_signal, manifold, rng)
        next_received_signal = propagate(received_signal, 1ms)
        measurement2 = @inferred get_measurement(4000, next_received_signal, manifold, rng)
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
        emitters = (sat1,)
        receiver = @inferred Receiver(
            4e6Hz,
            noise_std = 5.0,
            gain_phase_mism_crosstalk = SMatrix{2,2,ComplexF64}(I)
        )
        received_signal = @inferred ReceivedSignal(emitters, receiver)

        rng = MersenneTwister(1234)
        measurement = @inferred get_measurement(4000, received_signal, manifold, rng)
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


end
