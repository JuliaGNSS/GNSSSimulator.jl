@testset "Jammer" begin

    @testset "Create jammer with integers" begin
        jammer = CWJammer(
            1,
            12dB,
            doppler = 1100Hz,
            phase = 1,
            exists = true,
            doa = SVector(0, 0, 1)
        )
        @test @inferred(get_jammer_to_noise_ratio(jammer)) == 12dB
        @test @inferred(get_carrier_doppler(jammer)) == 1100Hz
    end

    @testset "CW Jammer" begin
        jammer = CWJammer(
            1,
            12dB,
            doppler = 1100.0Hz,
            phase = π,
            exists = true,
            doa = SVector(0, 0, 1)
        )
        @test @inferred(get_jammer_to_noise_ratio(jammer)) == 12dB
        @test @inferred(get_amplitude(jammer, 2/Hz, 1000Hz)) ≈ sqrt(10^(12/10) * 2000)
        @test @inferred(get_existence(jammer)) == true
        @test @inferred(get_id(jammer)) == 1
        @test @inferred(get_carrier_doppler(jammer)) == 1100Hz
        @test @inferred(get_carrier_phase(jammer)) ≈ π

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.gen_signal!(
            jammer,
            2.5e6Hz,
            100Hz,
            1/Hz,
            2500,
            rng
        )

        reference_signal = cis.(2π .* (0:2499) .* (1100Hz .+ 100Hz) ./ 2.5e6Hz .+ π) .*
            sqrt(10^(12/10) * 2.5e6)

        @test sqrt(sum(abs2.(signal .- reference_signal)) / 2500) < 1e-2 *
            sqrt(10^(12/10) * 2.5e6)

        next_jammer = @inferred GNSSSimulator.propagate(
            jammer,
            2500,
            100.0Hz,
            2.5e6Hz,
            Random.GLOBAL_RNG
        )
        @test @inferred(get_carrier_phase(next_jammer)) ≈ mod2pi(
            π + 2π * 1200.0Hz * 2500 / 2.5e6Hz + π
        ) - π
        @test @inferred(get_carrier_doppler(next_jammer)) == 1100.0Hz
        @test @inferred(get_jammer_to_noise_ratio(next_jammer)) == 12dB
        @test @inferred(get_existence(next_jammer)) == true
        @test @inferred(get_id(next_jammer)) == 1
        @test @inferred(get_doa(next_jammer)) == SVector(0.0, 0.0, 1.0)

        jammer = @inferred CWJammer(1, 10.0dB)
        @test get_id(jammer) == 1
        @test get_jammer_to_noise_ratio(jammer) == 10.0dB
    end

    @testset "Noise Jammer" begin
        noise_jammer = NoiseJammer(1, 12dB, exists = true, doa = SVector(0, 0, 1))
        @test @inferred(get_jammer_to_noise_ratio(noise_jammer)) == 12dB
        @test @inferred(get_amplitude(noise_jammer, 2/Hz, 1000Hz)) ≈
            sqrt(10^(12/10) * 2000)
        @test @inferred(get_existence(noise_jammer)) == true
        @test @inferred(get_id(noise_jammer)) == 1

        rng = MersenneTwister(1234)
        signal = @inferred GNSSSimulator.gen_signal!(
            noise_jammer,
            2.5e6Hz,
            120Hz,
            1/Hz,
            2500,
            rng
        )
        rng = MersenneTwister(1234)
        @test signal ≈ randn(rng, ComplexF64, 2500) .* sqrt(10^(12/10) * 2.5e6)

        next_jammer = @inferred GNSSSimulator.propagate(
            noise_jammer,
            2500,
            100.0Hz,
            2.5e6Hz,
            Random.GLOBAL_RNG
        )
        @test @inferred(get_jammer_to_noise_ratio(next_jammer)) == 12dB
        @test @inferred(get_existence(next_jammer)) == true
        @test @inferred(get_doa(next_jammer)) == SVector(0.0, 0.0, 1.0)
        @test @inferred(get_id(next_jammer)) == 1
    end
end
