const DEFAULT_DOAS = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
                     -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
                      0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]

srand(1234)
@testset "Existing sats" begin
    existing_sats = @inferred GNSSSimulator.sim_existing_sats([true, true, false])
    @test @inferred(existing_sats(1)) == [true, true, false]

    existing_sats_over_time = [true true true; true true false; true false false]
    gen_existing_sats = @inferred GNSSSimulator.sim_existing_sats(existing_sats_over_time, 10)

    @test @inferred(gen_existing_sats(0)) == [true, true, true]
    @test @inferred(gen_existing_sats(0.1)) == [true, true, false]
    @test @inferred(gen_existing_sats(0.2)) == [true, false, false]

end

@testset "Pseudo post corr signal" begin
    @testset "Power" begin
        pseudo_post_corr_signal = @inferred GNSSSimulator.sim_pseudo_post_corr_signal(1, 0)
        signal_ampl_and_phase = @inferred pseudo_post_corr_signal(0, true)
        @test abs(signal_ampl_and_phase) ≈ 1
    end

    @testset "Phase variance between satellites" begin
        gen_signal_ampl_and_phase = @inferred GNSSSimulator.sim_pseudo_post_corr_signal(1000, 0, 0.1)
        signal_ampl_and_phase = @inferred gen_signal_ampl_and_phase(0, trues(1000))
        @test var(angle.(signal_ampl_and_phase)) ≈ 0.1 atol = 0.01
    end
end

@testset "Noise" begin
    num_ants = 4
    num_sats = 100
    noise_power_dB = -15
    gen_noise = @inferred GNSSSimulator.sim_noise(noise_power_dB, num_ants)
    noise = @inferred(gen_noise(0, trues(num_sats)))
    @test noise * noise' / num_sats ≈ eye(num_ants) * 10^(noise_power_dB / 20) atol = 0.4
end

@testset "DOA over time" begin

    existing_sats = [true, true, true, true, true, true, true, true, true, true, false]
    gen_doas = @inferred GNSSSimulator.sim_doas()
    @test @inferred(gen_doas(1, existing_sats)) == DEFAULT_DOAS[:,1:10]

    num_sats = 2
    num_timestamps = 2
    doas_data = zeros(3,num_sats,num_timestamps)
    doas_data[:,1,1] = [1;0;0]
    doas_data[:,1,2] = [0;1;0]
    gen_doas = @inferred GNSSSimulator.sim_doas(doas_data, 10)

    @test @inferred(gen_doas(0.0, trues(2))) == [1 0;0 0;0 0]
    @test @inferred(gen_doas(0.1, trues(2))) == [0 0;1 0;0 0]

    @test @inferred(gen_doas(0.0, [true; false])) ≈ [1;0;0]
    @test @inferred(gen_doas(0.1, [true; false])) ≈ [0;1;0]
end

@testset "Attitude" begin

    attitude = @inferred GNSSSimulator.sim_attitude(0.0, 0.0, 0.0)
    @test @inferred(attitude(2)) == RotXYZ(0.0, 0.0, 0.0)

    attitude_over_time = [
            20 21; # Roll
            10 11; # Pitch
            120 121 # Yaw
        ] * π / 180
    attitude =  @inferred GNSSSimulator.sim_attitude(attitude_over_time, 10)
    @test @inferred(attitude(0)) ≈ RotXYZ(20 * π / 180, 10 * π / 180, 120 * π / 180)
    @test @inferred(attitude(0.1)) ≈ RotXYZ(21 * π / 180, 11 * π / 180, 121 * π / 180)

    attitude =  @inferred GNSSSimulator.sim_attitude([20; 10; 120] * π / 180, 1e-5, 1 * π / 180, 1 * π / 180, 1 * π / 180)
    rots = map(t -> RotXYZ(@inferred(attitude(t))), 1:1000)
    attitudes = hcat(map(T -> [T.theta1, T.theta2, T.theta3], rots)...)
    @test var(attitudes, 2) ≈ [1;1;1] * π / 180 atol = 0.5 * π / 180
end