@testset "Jammer" begin
    @testset "CW Jammer" begin
        jammer = CWJammer(1, 100.0Hz, 1.0, 10.0, true, SVector(0.0, 0.0, 1.0))

        next_jammer = propagate(jammer, 1μs)
        @test GNSSSimulator.get_carrier_phase(next_jammer) ≈ 1.0 + 2π * 100.0Hz * 1μs
        @test GNSSSimulator.get_carrier_doppler(next_jammer) == 100.0Hz
        @test GNSSSimulator.get_amplitude(next_jammer) == 10
        @test GNSSSimulator.get_existence(next_jammer) == true
        @test GNSSSimulator.get_doa(next_jammer) == SVector(0.0, 0.0, 1.0)

        get_steer_vec(doa, attitude) = 1 # No steering vector
        @test get_signal(jammer, nothing, get_steer_vec) ≈ cis(1.0) * 10
    end

    @testset "Noise Jammer" begin
        noise = randn(ComplexF64)
        jammer = NoiseJammer(1, noise, 10.0, true, SVector(0.0, 0.0, 1.0))

        next_jammer = propagate(jammer, 1μs)
        @test GNSSSimulator.get_noise(jammer) == noise
        @test GNSSSimulator.get_amplitude(next_jammer) == 10
        @test GNSSSimulator.get_existence(next_jammer) == true
        @test GNSSSimulator.get_doa(next_jammer) == SVector(0.0, 0.0, 1.0)

        get_steer_vec(doa, attitude) = 1 # No steering vector
        @test get_signal(jammer, nothing, get_steer_vec) == noise * 10
    end
end
