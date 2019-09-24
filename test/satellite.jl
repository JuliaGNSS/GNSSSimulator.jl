@testset "Satellite" begin
    @testset "Auxiliarly functions" begin
        doa = Spherical(1.0, 0.0, π / 2) # Zenith
        distance_from_earth_center = EARTH_RADIUS + 40000m
        distance = @inferred GNSSSimulator.calc_sat_user_distance(distance_from_earth_center, doa)

        @test distance == 40000m

        @test @inferred(GNSSSimulator.calc_carrier_phase(0m, 20Hz)) == 0

        @test @inferred(GNSSSimulator.calc_doppler(distance_from_earth_center, doa, 14_000.0m / 3.6s, 1.157542e6Hz)) ≈ 0Hz rtol = 1
    end

    @testset "Satellite signal" begin
        sat = ConstantDopplerSatellite(GPSL1, 1, carrier_doppler = 1000.0Hz, carrier_phase = π / 2, code_phase = 100.0, cn0 = 45dBHz)
        next_sat = @inferred propagate(sat, 1/2e6Hz)
        @test @inferred(GNSSSimulator.get_carrier_doppler(next_sat)) == 1000Hz
        @test @inferred(GNSSSimulator.get_code_doppler(next_sat)) ≈ 1000Hz / 1540
        @test @inferred(GNSSSimulator.get_carrier_phase(next_sat)) ≈ π / 2 + 2π * 1000Hz / 2e6Hz
        @test @inferred(GNSSSimulator.get_code_phase(next_sat)) ≈ 100 + (1023e3Hz + 1000Hz / 1540) / 2e6Hz
        @test @inferred(GNSSSimulator.get_existence(next_sat)) == true
        @test @inferred(GNSSSimulator.get_amplitude(next_sat)) == 10^(45 / 20)
        @test @inferred(GNSSSimulator.get_prn(next_sat)) == 1

        manifold = IdealManifold()

        signal = @inferred get_signal(next_sat, nothing, manifold)
        @test signal ≈ cis(π / 2 + 2π * 1000Hz / 2e6Hz) * get_code(GPSL1, 100 + (1023e3Hz + 1000Hz / 1540) / 2e6Hz, 1) * 10^(45 / 20)
    end
end
