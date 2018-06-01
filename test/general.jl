srand(1234)
@testset "Existing sats" begin
    existing_sats_over_time = [true true true; true true false; true false false]
    gen_existing_sats = @inferred GNSSSimulator.init_gen_existing_sats(existing_sats_over_time, 10)

    @test @inferred gen_existing_sats(0) == [true, true, true]
    @test @inferred gen_existing_sats(0.1) == [true, true, false]
    @test @inferred gen_existing_sats(0.2) == [true, false, false]

end

@testset "Signal ampl and phase" begin
    @testset "Power" begin
        gen_signal_ampl_and_phase = @inferred GNSSSimulator.init_gen_signal_ampl_and_phase(1, 0, 0, 0, 0)
        signal_ampl_and_phase = @inferred gen_signal_ampl_and_phase(0, true)
        @test abs(signal_ampl_and_phase) ≈ 1
    end

    @testset "Amplitude and phase standard deviation between satellites" begin
        gen_signal_ampl_and_phase = @inferred GNSSSimulator.init_gen_signal_ampl_and_phase(1000, 0.1, 0.1, 0, 0)
        signal_ampl_and_phase = @inferred gen_signal_ampl_and_phase(0, trues(1000))
        @test std(abs.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
        @test std(angle.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
    end

    @testset "Amplitude and phase standard deviation over time" begin
        gen_signal_ampl_and_phase = @inferred GNSSSimulator.init_gen_signal_ampl_and_phase(1, 0, 0, 0.1, 0.1)
        @inferred gen_signal_ampl_and_phase(1, true)
        signal_ampl_and_phase = map(t -> gen_signal_ampl_and_phase(t, true), 1:1000)
        @test std(abs.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
        @test std(angle.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
    end

end

@testset "Noise" begin
    num_ants = 4
    num_sats = 100
    noise_power_dB = -15
    gen_gen_noise = @inferred GNSSSimulator.init_gen_noise(noise_power_dB, num_ants)
    noise = @inferred gen_gen_noise(0, trues(num_sats))
    @test noise * noise' / num_sats ≈ eye(num_ants) * 10^(noise_power_dB / 20) atol = 0.4
end
