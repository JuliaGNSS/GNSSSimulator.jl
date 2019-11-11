@testset "Receiver" begin

    receiver = Receiver(4e6Hz, 0.0Hz, SMatrix{4,4}(I), RotXYZ(0.0, 0.0, 0.0), 0.0)
    next_receiver = GNSSSimulator.propagate(receiver, 1ms, Random.GLOBAL_RNG)
    @test get_gain_phase_mism_crosstalk(next_receiver) == SMatrix{4,4}(I)
    @test get_attitude(next_receiver) == RotXYZ(0.0, 0.0, 0.0)
    @test get_noise_std(next_receiver) == 0.0
    @test get_sample_frequency(next_receiver) == 4e6Hz
    @test get_intermediate_frequency(next_receiver) == 0.0Hz
    @test get_gain_phase_mism_crosstalk(next_receiver) == SMatrix{4,4}(I)

    receiver = Receiver(
        4e6Hz,
        intermediate_frequency = 0.0Hz,
        gain_phase_mism_crosstalk = SMatrix{4,4,ComplexF64}(I),
        attitude = RotXYZ(0.0, 0.0, 0.0),
        n0 = 1/Hz
    )
    @test get_gain_phase_mism_crosstalk(receiver) == SMatrix{4,4}(I)
    @test get_attitude(receiver) == RotXYZ(0.0, 0.0, 0.0)
    @test get_noise_std(receiver) == sqrt(4e6)
    @test get_sample_frequency(receiver) == 4e6Hz
    @test get_intermediate_frequency(receiver) == 0.0Hz
    @test get_gain_phase_mism_crosstalk(receiver) == SMatrix{4,4}(I)
end
