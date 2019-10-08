@testset "Existence" begin
    @testset "Static Existence" begin
        @test GNSSSimulator.propagate(true, 1s, Random.GLOBAL_RNG) == true
        @test GNSSSimulator.propagate(false, 1s, Random.GLOBAL_RNG) == false
        @test get_existence(true) == true
        @test get_existence(false) == false
    end
    @testset "Dynamic Existence" begin
        existences = [true, false, true, false, true, false]
        dynamic_existence = @inferred DynamicExistence(existences, 1.0Hz, 0.0s)

        next_dynamic_existence = @inferred GNSSSimulator.propagate(dynamic_existence, 1s, Random.GLOBAL_RNG)
        @test next_dynamic_existence == DynamicExistence(existences, 1.0Hz, 1.0s)
        @test get_existence(next_dynamic_existence) == false
        @test get_existence(GNSSSimulator.propagate(next_dynamic_existence, 3s, Random.GLOBAL_RNG)) == true

    end
end
