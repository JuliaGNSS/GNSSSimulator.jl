@testset "Gain and phase mismatch and crosstalk" begin

    @testset "Crosstalk amplitude" begin
        crosstalk_dB = -10
        gen_gp_mism_crosstalk = @inferred GNSSSimulator.sim_gain_phase_mism_and_crosstalk(4, crosstalk_dB, 0, 0, 0, 0)

        gp_mism_crosstalk = @inferred gen_gp_mism_crosstalk(0)

        @test abs(gp_mism_crosstalk[1,1]) == abs(gp_mism_crosstalk[2,2])
        @test abs(gp_mism_crosstalk[1,1]) * 10^(crosstalk_dB / 10) ≈ abs(gp_mism_crosstalk[1,2])
        @test abs(gp_mism_crosstalk[1,2]) == abs(gp_mism_crosstalk[2,1])
    end

    @testset "Normalization" begin
        normed_matrix = @inferred GNSSSimulator.normalize_gain_phase_mism_and_crosstalk(ones(4,4))
        @test normed_matrix == ones(4,4) / 2
    end
end

@testset "Steering vector" begin

    doas = [1 0; 0 0; 0 1]
    attitude = RotXYZ(0.0, 0.0, 1.0 * π)
    gen_steering_vectors = @inferred GNSSSimulator.sim_steering_vectors(a -> [a[1] + 0.0im, a[1] + 0.0im, a[2] + 0.0im, a[3] + 0.0im])
    @test @inferred(gen_steering_vectors(0, attitude, doas)) ≈ [-1.0 0.0; -1.0 0.0; 0.0 0.0; 0.0 1.0]
end

@testset "Measurement" begin

    doas = sim_doas()
    existing_sats = sim_existing_sats(trues(11))
    pseudo_post_corr_signal = sim_pseudo_post_corr_signal(11, 0)
    attitude = sim_attitude(0.0, 0.0, 0.0)
    gain_phase_mism_and_crosstalk = sim_gain_phase_mism_and_crosstalk(4, -15)
    steering_vectors = sim_steering_vectors(a -> [a[1] + 0.0im, a[1] + 0.0im, a[2] + 0.0im, a[3] + 0.0im])
    noise = sim_noise(-15, 4)

    measurement = @inferred sim_post_corr_measurement(
        existing_sats,
        pseudo_post_corr_signal,
        attitude,
        doas,
        gain_phase_mism_and_crosstalk,
        steering_vectors,
        noise)

    𝐘, internal_states = @inferred measurement(0)
    @test size(𝐘) == (4,11)
end
