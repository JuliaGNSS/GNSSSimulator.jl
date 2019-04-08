@testset "System" begin

    enu_doa = Spherical(1.0, 0.0, π / 2) # Zenith
    distance_from_earth_center = EARTH_RADIUS + 40000m
    distance = calc_init_sat_user_distance(distance_from_earth_center, enu_doa)

    @test distance == 40000m

    @test GNSSSimulator.calc_carrier_phase(0m, 20Hz) == 0

    @test GNSSSimulator.calc_init_doppler(distance_from_earth_center, enu_doa, 14_000.0m / 3.6s, 1.157542e6Hz) ≈ 0Hz rtol = 1

    @test GNSSSimulator.doppler(0m / 1s, 20Hz) ≈ 0Hz
end
