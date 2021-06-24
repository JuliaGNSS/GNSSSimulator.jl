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

    @testset "Create satellite with integers" begin
        system = GPSL1()
        sat = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000Hz,
            carrier_phase = 1,
            code_phase = 100,
            cn0 = 45dBHz
        )

        @test @inferred(GNSSSimulator.get_carrier_phase(sat)) == 1
        @test @inferred(GNSSSimulator.get_code_phase(sat)) == 100
        @test @inferred(GNSSSimulator.get_carrier_doppler(sat)) == 1000Hz
    end

    @testset "Satellite signal" begin
        system = GPSL1()
        sat = ConstantDopplerSatellite(
            system,
            1,
            carrier_doppler = 1000.0Hz,
            carrier_phase = π / 2,
            code_phase = 100.0,
            cn0 = 45dBHz
        )

        @test @inferred(GNSSSimulator.get_carrier_phase(sat)) == π / 2
        @test @inferred(GNSSSimulator.get_code_phase(sat)) == 100

        code = Vector{Int8}(undef, 2500)

        code = @inferred GNSSSimulator.gen_code!(
            code,
            system,
            1023e3Hz + 1Hz,
            2.5e6Hz,
            120,
            1
        )

        @test code == get_code.(
            system,
            (0:2499) .* (1023e3 .+ 1) ./ 2.5e6 .+ 120,
            1
        )

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.gen_signal!(
            sat,
            2.5e6Hz,
            10.0Hz,
            1/Hz,
            2500,
            rng
        )

        code_doppler = 1000.0 * 1023e3 / 1.57542e9
        reference_signal = get_code.(
            system,
            (0:2499) .* (1023e3 .+ code_doppler) ./ 2.5e6 .+ 100,
            1
        ) .* cis.(2π .* (0:2499) .* (1000.0 + 10.0) ./ 2.5e6 .+ π / 2) .* sqrt(10^(45 / 10))
        @test sqrt(sum(abs2.(signal .- reference_signal)) / 2500) < 1e-2 * sqrt(10^(45 / 10))

        next_sat = @inferred GNSSSimulator.propagate(
            sat,
            2500,
            100.0Hz,
            2.5e6Hz,
            1.0/Hz,
            Random.GLOBAL_RNG
        )

        @test @inferred(get_carrier_doppler(next_sat)) == 1000Hz
        @test @inferred(get_code_doppler(next_sat)) ≈ 1000Hz / 1540
        @test @inferred(get_carrier_phase(next_sat)) ≈ mod2pi(
            π / 2 + 2π * 1100.0Hz * 2500 / 2.5e6Hz + π
        ) - π
        @test @inferred(get_code_phase(next_sat)) ≈ mod(
            100.0 + (1023e3Hz + 1000.0Hz / 1540) * 2500 / 2.5e6Hz,
            1023
        )
        @test @inferred(get_existence(next_sat)) == true
        @test @inferred(get_carrier_to_noise_density_ratio(next_sat)) == 45dBHz
        @test @inferred(get_prn(next_sat)) == 1
        @test @inferred(get_gnss_system(next_sat)) == system

    end
end
