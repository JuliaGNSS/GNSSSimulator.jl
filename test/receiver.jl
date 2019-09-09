@testset "Receiver" begin

    receiver = Receiver(SMatrix{4,4}(I), RotXYZ(0.0, 0.0, 0.0), 0.0)
    @test propagate(receiver, 1ms) == Receiver(SMatrix{4,4}(I), RotXYZ(0.0, 0.0, 0.0), 0.0)

    @test get_gain_phase_mism_crosstalk(receiver) == SMatrix{4,4}(I)
end
