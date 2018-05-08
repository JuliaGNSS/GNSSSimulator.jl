srand(1234)
@testset "Existing sats" begin

    existing_sats_over_time = [true true true; true true false; true false false]

    gen_existing_sats = GNSSSimulator.init_gen_existing_sats(existing_sats_over_time, 10)

    @test gen_existing_sats(0) == [true, true, true]
    @test gen_existing_sats(0.1) == [true, true, false]
    @test gen_existing_sats(0.2) == [true, false, false]

end

@testset "Signal ampl and phase" begin

    @testset "SNR" begin
        SNR_dB = 15
        gen_signal_ampl_and_phase = GNSSSimulator.init_gen_signal_ampl_and_phase(1, SNR_dB, 0, 0, 0, 0)
        signal_ampl_and_phase = gen_signal_ampl_and_phase(0, true)
        @test abs(signal_ampl_and_phase) == 10^(SNR_dB / 20)
    end

    @testset "Amplitude and phase standard deviation between satellites" begin
        gen_signal_ampl_and_phase = GNSSSimulator.init_gen_signal_ampl_and_phase(1000, 0, 0.1, 0.1, 0, 0)
        signal_ampl_and_phase = gen_signal_ampl_and_phase(0, trues(1000))
        @test std(abs.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
        @test std(angle.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
    end

    @testset "Amplitude and phase standard deviation over time" begin
        gen_signal_ampl_and_phase = GNSSSimulator.init_gen_signal_ampl_and_phase(1, 0, 0, 0, 0.1, 0.1)
        signal_ampl_and_phase = map(t -> gen_signal_ampl_and_phase(t, true), 1:1000)
        @test std(abs.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
        @test std(angle.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
    end

end

@testset "Noise" begin

    num_ants = 4
    num_sats = 100
    gen_gen_noise = GNSSSimulator.init_gen_noise(num_ants)
    noise = gen_gen_noise(0, trues(num_sats))
    @test noise * noise' / num_sats ≈ eye(num_ants) atol = 0.4

end
