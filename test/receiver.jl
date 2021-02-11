@testset "Receiver" begin

    receiver = Receiver(
        4e6Hz,
        intermediate_frequency = 0.0Hz,
        gain_phase_mism_crosstalk = SMatrix{4,4,ComplexF64}(I),
        attitude = RotXYZ(0.0, 0.0, 0.0),
        n0 = 1/Hz
    )
    next_receiver = GNSSSimulator.propagate(receiver, 1000, Random.GLOBAL_RNG)
    @test get_gain_phase_mism_crosstalk(next_receiver) == SMatrix{4,4}(I)
    @test get_attitude(next_receiver) == RotXYZ(0.0, 0.0, 0.0)
    @test get_noise_density(next_receiver) == 1/Hz
    @test get_sampling_frequency(next_receiver) == 4e6Hz
    @test get_intermediate_frequency(next_receiver) == 0.0Hz
    @test get_gain_phase_mism_crosstalk(next_receiver) == SMatrix{4,4}(I)
    @test @inferred(get_noise_std(next_receiver)) == sqrt(1/Hz * 4e6Hz)

    receiver = Receiver(
        4e6Hz,
        intermediate_frequency = 0.0Hz,
        gain_phase_mism_crosstalk = SMatrix{4,4,ComplexF64}(I),
        attitude = RotXYZ(0.0, 0.0, 0.0),
        n0 = 1/Hz
    )
    @test get_gain_phase_mism_crosstalk(receiver) == SMatrix{4,4}(I)
    @test get_attitude(receiver) == RotXYZ(0.0, 0.0, 0.0)
    @test get_noise_density(receiver) == 1/Hz
    @test get_sampling_frequency(receiver) == 4e6Hz
    @test get_intermediate_frequency(receiver) == 0.0Hz
    @test get_gain_phase_mism_crosstalk(receiver) == SMatrix{4,4}(I)
end
