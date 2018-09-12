@testset "Existence" begin
    static_exist = true 
    dyn_exist = GNSSSimulator.DynamicExistence([false; true; false; true], 1Hz)

    @test true == sim_existence(3s, static_exist)
    @test true == sim_existence(1s, dyn_exist)
    @test true == sim_existence(5s, dyn_exist) # time exceeds data length --> return last available value
end