const CART_COORD = SVector{3}([1.0; 2.0; 3.0])
const SPH_COORD = SphericalFromCartesian()(CART_COORD)
const CART_COORD_VEC = [[CART_COORD]; [2 * CART_COORD]; [3 * CART_COORD]]
const SPH_COORD_VEC = SphericalFromCartesian().(CART_COORD_VEC)
const DYN_DOA_CART = GNSSSimulator.DynamicDOA(CART_COORD_VEC, 1Hz)
const DYN_DOA_SPH = GNSSSimulator.DynamicDOA(SPH_COORD_VEC, 1Hz) 

@testset "DOAs" begin
    @testset "Static DOA" begin
        @test CART_COORD == sim_doa(0s, SVector{3}(CART_COORD))
        @test CART_COORD == sim_doa(0s, SPH_COORD)
    end
    
   @testset "Dynamic DOA" begin
        @test CART_COORD_VEC[2] == sim_doa(1s, DYN_DOA_CART) 
        @test CART_COORD_VEC[3] == sim_doa(5s, DYN_DOA_CART)  # time exceeds data length --> return last available value

        @test CART_COORD_VEC[2] == sim_doa(1s, DYN_DOA_SPH)
        @test CART_COORD_VEC[3] ≈ sim_doa(5s, DYN_DOA_SPH) # time exceeds data length --> return last available value
    end

    @testset "Linear Dynamic DOA"  begin
        lin_dyn_doa = GNSSSimulator.LinearDynamicDOA(SPH_COORD, 0.1Hz, 0.2Hz)
        lin_dyn_doa_with_default = GNSSSimulator.LinearDynamicDOA(init_DOA = SPH_COORD)

        @test CartesianFromSpherical()(Spherical(SPH_COORD.r, SPH_COORD.θ + 0.1 * 3, SPH_COORD.ϕ + 0.2 * 3)) == sim_doa(3s, lin_dyn_doa)
        @test CART_COORD_VEC[1] == sim_doa(1s, lin_dyn_doa_with_default)
    end
end