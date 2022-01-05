@testset "Synthetic satellite" begin
    system = GPSL1()
    sat = ConstantDopplerSatellite(
        system,
        1,
        carrier_doppler = 1000.0Hz,
        carrier_phase = π / 2,
        code_phase = 100.0,
        cn0 = 45dBHz
    )
    synthetic_sat = SyntheticSatellite(sat)

    @test @inferred(get_carrier_phase(synthetic_sat)) == π / 2
    @test @inferred(get_code_phase(synthetic_sat)) == 100.0

    rng = MersenneTwister(1234)
    signal = @inferred GNSSSimulator.gen_signal!(
        synthetic_sat,
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
    ) .* cis.(2π .* (0:2499) .* (1000.0 + 10.0) ./ 2.5e6 .+ π / 2) .*
        sqrt(10^(45 / 10))
    @test sqrt(sum(abs2.(signal .- reference_signal)) / 2500) < 1e-2 * sqrt(10^(45 / 10))

    next_synthetic_sat = @inferred GNSSSimulator.propagate(
        synthetic_sat,
        1,
        100.0Hz,
        2e6Hz,
        1.0/Hz,
        Random.GLOBAL_RNG
    )
    @test @inferred(get_carrier_doppler(next_synthetic_sat)) == 1000Hz
    @test @inferred(get_code_doppler(next_synthetic_sat)) ≈ 1000Hz / 1540
    @test @inferred(get_carrier_phase(next_synthetic_sat)) ≈ π / 2 + 2π * 1100Hz / 2e6Hz
    @test @inferred(get_code_phase(next_synthetic_sat)) ≈ 100 + (1023e3Hz + 1000Hz / 1540) / 2e6Hz
    @test @inferred(get_existence(next_synthetic_sat)) == true
    @test @inferred(get_carrier_to_noise_density_ratio(next_synthetic_sat)) ≈ 45dBHz
    @test @inferred(get_prn(next_synthetic_sat)) == 1
    @test @inferred(get_gnss_system(next_synthetic_sat)) == system

    manifold = IdealManifold(
        1575420e3,
        0.1904 / 4 * SVector(SVector(1, 1, 0), SVector(-1, 1, 0)),
    )
    @test get_steer_vec(manifold, synthetic_sat, RotXYZ(0.0,0.0,0.0)) == [1.0, 1.0]

    synthetic_sat_with_steer_vec = SyntheticSatellite(sat, SVector(1.0, 1.0im))
    @test get_steer_vec(manifold, synthetic_sat_with_steer_vec, RotXYZ(0.0,0.0,0.0)) == [1.0, 1.0im]

end
