@testset "Receiver" begin

    receiver = Receiver(SMatrix{4,4}(I), RotXYZ(0.0, 0.0, 0.0), 0.0)
    next_receiver = propagate(receiver, 1ms)
    @test get_gain_phase_mism_crosstalk(next_receiver) == SMatrix{4,4}(I)
    @test get_attitude(next_receiver) == RotXYZ(0.0, 0.0, 0.0)
    @test get_noise(next_receiver) ≈ SVector{4, ComplexF64}(0.0, 0.0, 0.0, 0.0)

    @test get_gain_phase_mism_crosstalk(receiver) == SMatrix{4,4}(I)
end
