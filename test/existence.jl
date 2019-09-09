@testset "Existence" begin
    @testset "Static Existence" begin
        @test propagate(true, 1s) == true
        @test propagate(false, 1s) == false
        @test get_existence(true) == true
        @test get_existence(false) == false
    end
    @testset "Dynamic Existence" begin
        existences = [true, false, true, false, true, false]
        dynamic_existence = DynamicExistence(existences, 0.0s, 1.0Hz)

        next_dynamic_existence = propagate(dynamic_existence, 1s)
        @test next_dynamic_existence == DynamicExistence(existences, 1.0s, 1.0Hz)
        @test get_existence(next_dynamic_existence) == false
        @test get_existence(propagate(next_dynamic_existence, 3s)) == true

    end
end
