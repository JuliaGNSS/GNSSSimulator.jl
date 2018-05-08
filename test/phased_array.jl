@testset "Gain and phase mismatch and crosstalk" begin

    # Tests on standard deviation are missing

    @testset "Crosstalk amplitude" begin
        crosstalk_dB = -10
        gen_gp_mism_crosstalk = GNSSSimulator.init_gen_gain_and_phase_mism_and_crosstalk(2, 0, 0, 0, 0, crosstalk_dB, 0, 0, 0, 0)

        gp_mism_crosstalk = gen_gp_mism_crosstalk(0)

        @test abs(gp_mism_crosstalk[1,1]) == abs(gp_mism_crosstalk[2,2])
        @test abs(gp_mism_crosstalk[1,1]) * 10^(crosstalk_dB / 10) ≈ abs(gp_mism_crosstalk[1,2])
        @test abs(gp_mism_crosstalk[1,2]) == abs(gp_mism_crosstalk[2,1])
    end

    @testset "Gain" begin
        gen_gp_mism_crosstalk = GNSSSimulator.init_gen_gain_and_phase_mism_and_crosstalk(2, 0, 0, 0, 0, -Inf, 0, 0, 0, 0)

        gp_mism_crosstalk = gen_gp_mism_crosstalk(0)

        @test abs(gp_mism_crosstalk[1,1]) == abs(gp_mism_crosstalk[2,2])
        @test abs(gp_mism_crosstalk[1,1]) == 1
    end

    @testset "Normalization" begin
        normed_matrix = GNSSSimulator.normalize_gain_and_phase_mism_and_crosstalk(ones(4,4))
        @test mapslices(norm, normed_matrix, 1) == [1.0 1.0 1.0 1.0]
    end
end

@testset "DOA over time" begin
    num_sats = 2
    num_timestamps = 2
    doas = zeros(3,num_sats,num_timestamps)
    doas[:,1,1] = [1;0;0]
    doas[:,1,2] = [0;1;0]
    gen_doas = GNSSSimulator.init_gen_doas(doas, 10)

    @test gen_doas(0.0, trues(2)) == [1 0;0 0;0 0]
    @test gen_doas(0.1, trues(2)) == [0 0;1 0;0 0]

    @test gen_doas(0.0, [true; false]) ≈ [1;0;0]
    @test gen_doas(0.1, [true; false]) ≈ [0;1;0]
end

@testset "Attitude" begin

    @testset "Attitude" begin
        attitude_over_time = [
            20 21; # Roll
            10 11; # Pitch
            120 121 # Yaw
        ] * π / 180
        gen_attitude = GNSSSimulator.init_gen_attitude(attitude_over_time, 10, 0)
        @test gen_attitude(0) == RotXYZ(20 * π / 180, 10 * π / 180, 120 * π / 180)
        @test gen_attitude(0.1) == RotXYZ(21 * π / 180, 11 * π / 180, 121 * π / 180)
    end

    @testset "Attitude over time" begin
        attitude = [20; 10; 120] * π / 180
        gen_attitude = GNSSSimulator.init_gen_attitude(attitude, 1/10000, 1 * π / 180)
        rots = map(t -> RotXYZ(gen_attitude(t)), 1:1000)
        attitudes = hcat(map(T -> [T.theta1, T.theta2, T.theta3], rots)...)
        @test std(attitudes, 2) ≈ [1;1;1] * π / 180 atol = 0.5 * π / 180
    end
end

@testset "Steering vector" begin

    doas = [1 0; 0 0; 0 1]
    attitude = RotXYZ(0.0, 0.0, 1.0 * π)
    gen_steering_vectors = GNSSSimulator.init_gen_steering_vectors(a -> [a.ϕ;a.ϕ;a.θ;a.θ])
    @test gen_steering_vectors(0, attitude, doas) ≈ [0.0 π / 2; 0.0 π / 2; π 0.0; π 0.0]
end
