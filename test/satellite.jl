@testset "Satellite" begin
    @testset "Auxiliarly functions" begin
        doa = Spherical(1.0, 0.0, π / 2) # Zenith
        distance_from_earth_center = EARTH_RADIUS + 40000m
        distance = @inferred GNSSSimulator.calc_sat_user_distance(
            distance_from_earth_center,
            doa
        )

        @test distance == 40000m

        @test @inferred(GNSSSimulator.calc_carrier_phase(0m, 20Hz)) == 0

        doppler = @inferred GNSSSimulator.calc_doppler(
            distance_from_earth_center,
            doa,
            14_000.0m / 3.6s,
            1.157542e6Hz
        )
        @test doppler ≈ 0Hz rtol = 1
    end

    @testset "Phase wrap" begin
        sat = ConstantDopplerSatellite(GPSL1, 1)
        phase = GNSSSimulator.SatellitePhase(0.25, 1023.25)
        phase_wrap = GNSSSimulator.SatellitePhaseWrap(0, 0)

        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(phase_wrap, phase, sat)
        @test next_phase_wrap == GNSSSimulator.SatellitePhaseWrap(0, 1023)

        phase = GNSSSimulator.SatellitePhase(0.55, 1023.55)
        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(
            next_phase_wrap,
            phase,
            sat
        )
        @test next_phase_wrap == GNSSSimulator.SatellitePhaseWrap(1, 1023)
    end

    @testset "Satellite signal" begin
        sat = ConstantDopplerSatellite(
            GPSL1,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )

        @test @inferred(GNSSSimulator.get_carrier_phase(sat)) == π / 2
        @test @inferred(GNSSSimulator.get_code_phase(sat)) == 100

        phase = @inferred GNSSSimulator.calc_phase(
            sat,
            0,
            100.0Hz,
            2e6Hz
        )
        @test phase == GNSSSimulator.SatellitePhase(0.25, 100.0)

        phase_wrap = @inferred GNSSSimulator.init_phase_wrap(sat)
        @test phase_wrap == GNSSSimulator.SatellitePhaseWrap(0, 0)

        signal = @inferred GNSSSimulator.get_signal(
            sat,
            phase,
            phase_wrap,
            1.0 + 0.0im,
            Random.GLOBAL_RNG
        )
        @test signal ≈ 1.0 * 10^(45 / 20) * cis(π / 2) * get_code(GPSL1, 100.0, 1)
        # cis(π / 2) == cis_vfast(π / 2)

        signal = @inferred GNSSSimulator.get_signal(
            sat,
            phase,
            phase_wrap,
            SVector(1.0im, 2.0im),
            Random.GLOBAL_RNG
        )
        @test signal ≈ SVector(1im, 2im) * 10^(45 / 20) * cis(π / 2) *
            get_code(GPSL1, 100.0, 1) # cis(π / 2) == cis_vfast(π / 2)

        next_phase = @inferred GNSSSimulator.calc_phase(
            sat,
            1,
            100.0Hz,
            2e6Hz
        )
        @test next_phase.carrier == 0.25 + 1100.0Hz / 2e6Hz
        @test next_phase.code == 100.0 + (1023e3Hz + 1000.0Hz / 1540) / 2e6Hz

        next_phase_wrap = @inferred GNSSSimulator.update_phase_wrap(
            phase_wrap,
            next_phase,
            sat
        )

        @test next_phase_wrap == GNSSSimulator.SatellitePhaseWrap(0, 0)

        signal = @inferred GNSSSimulator.get_signal(
            sat,
            next_phase,
            next_phase_wrap,
            1.0 + 0.0im,
            Random.GLOBAL_RNG
        )
        @test signal ≈ 1.0 * 10^(45 / 20) *
            GNSSSignals.cis_vfast(π / 2 + 2π * 1100.0Hz / 2e6Hz) *
            get_code(GPSL1, 100.0 + (1023e3Hz + 1000.0Hz / 1540) / 2e6Hz, 1)

        next_sat = @inferred GNSSSimulator.propagate(
            sat,
            100.0Hz,
            1/2e6Hz,
            Random.GLOBAL_RNG
        )

        @test @inferred(get_carrier_doppler(next_sat)) == 1000Hz
        @test @inferred(get_code_doppler(next_sat)) ≈ 1000Hz / 1540
        @test @inferred(get_carrier_phase(next_sat)) ≈ π / 2 + 2π * 1100.0Hz / 2e6Hz
        @test @inferred(get_code_phase(next_sat)) ≈ 100.0 + (1023e3Hz + 1000.0Hz / 1540) /
            2e6Hz
        @test @inferred(get_existence(next_sat)) == true
        @test @inferred(get_amplitude(next_sat)) == 10^(45 / 20)
        @test @inferred(get_prn(next_sat)) == 1
        @test @inferred(get_gnss_system(next_sat)) == GPSL1


    end
end
