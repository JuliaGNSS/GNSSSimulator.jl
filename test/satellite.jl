@testset "Satellite" begin
    @testset "Auxiliarly functions" begin
        doa = Spherical(1.0, 0.0, π / 2) # Zenith
        distance_from_earth_center = EARTH_RADIUS + 40000m
        distance = GNSSSimulator.calc_sat_user_distance(distance_from_earth_center, doa)

        @test distance == 40000m

        @test GNSSSimulator.calc_carrier_phase(0m, 20Hz) == 0

        @test GNSSSimulator.calc_doppler(distance_from_earth_center, doa, 14_000.0m / 3.6s, 1.157542e6Hz) ≈ 0Hz rtol = 1
    end

    @testset "Satellite signal" begin
        gpsl1 = GPSL1()
        sat = ConstantDopplerSatellite(1, gpsl1, carrier_doppler = 1000.0Hz, carrier_phase = π / 2, code_phase = 100.0, cn0 = 45dBHz)
        next_sat = propagate(sat, 1/2e6Hz)
        @test GNSSSimulator.get_carrier_doppler(next_sat) == 1000Hz
        @test GNSSSimulator.get_code_doppler(next_sat) ≈ 1000Hz / 1540
        @test GNSSSimulator.get_carrier_phase(next_sat) ≈ π / 2 + 2π * 1000Hz / 2e6Hz
        @test GNSSSimulator.get_code_phase(next_sat) ≈ 100 + 1000Hz / 1540 / 2e6Hz
        @test GNSSSimulator.get_existence(next_sat) == true
        @test GNSSSimulator.get_amplitude(next_sat) == 10^(45 / 20)
        @test GNSSSimulator.get_prn(next_sat) == 1
        @test GNSSSimulator.get_system(next_sat) == gpsl1

        get_steer_vec(doa, attitude) = 1 # No steering vector

        signal = get_signal(next_sat, nothing, get_steer_vec)
        @test signal ≈ cis(π / 2 + 2π * 1000Hz / 2e6Hz) * gpsl1.codes[1 + floor(Int, 100 + 1000Hz / 1540 / 2e6Hz), 1] * 10^(45 / 20)
    end
end
