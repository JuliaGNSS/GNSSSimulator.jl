@testset "Structural interference" begin
    system = GPSL1()

    sat = ConstantDopplerSatellite(
        system,
        1,
        carrier_doppler = 1000.0Hz,
        carrier_phase = π / 2,
        code_phase = 100.0,
        cn0 = 45dBHz
    )
    si = ConstantDopplerStructuralInterference(
        sat,
        -3dB,
        added_carrier_doppler = 10Hz,
        added_carrier_phase = π / 8,
        added_code_phase = 10.0
    )

    @test @inferred(get_carrier_phase(si)) == π / 2 + π / 8
    @test @inferred(get_code_phase(si)) == 110.0

    rng = MersenneTwister(1234)
    signal = @inferred GNSSSimulator.gen_signal!(
        si,
        2.5e6Hz,
        10.0Hz,
        1/Hz,
        2500,
        rng
    )

    code_doppler = 1000.0 * 1023e3 / 1.57542e9
    reference_signal = get_code.(
        system,
        (0:2499) .* (1023e3 .+ code_doppler) ./ 2.5e6 .+ 100 .+ 10.0,
        1
    ) .* cis.(2π .* (0:2499) .* (1000.0 + 10.0 + 10.0) ./ 2.5e6 .+ π / 2 .+ π / 8) .*
        sqrt(10^(42 / 10))
    @test sqrt(sum(abs2.(signal .- reference_signal)) / 2500) < 1e-2 * sqrt(10^(42 / 10))

    next_si = @inferred GNSSSimulator.propagate(si, 1, 100.0Hz, 2e6Hz, Random.GLOBAL_RNG)
    @test @inferred(get_carrier_doppler(next_si)) == 1010Hz
    @test @inferred(get_code_doppler(next_si)) ≈ 1010Hz / 1540
    @test @inferred(get_carrier_phase(next_si)) ≈ π / 2 + π / 8 + 2π * 1110Hz / 2e6Hz
    @test @inferred(get_code_phase(next_si)) ≈ 110 + (1023e3Hz + 1010Hz / 1540) / 2e6Hz
    @test @inferred(get_existence(next_si)) == true
    @test @inferred(get_carrier_to_noise_density_ratio(next_si)) ≈ 42dBHz
    @test @inferred(get_prn(next_si)) == 1
    @test @inferred(get_gnss_system(next_si)) == system
end
