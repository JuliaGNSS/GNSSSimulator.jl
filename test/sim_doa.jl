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

    @testset "Satellite channels with static spherical DOAs" begin
        channel_stat_signal_doa_sph = GNSSSimulator.SatelliteChannel(1, SPH_COORD, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, SVector{3}(CART_COORD), 0.0 + 0.0im, false)
        channel_stat_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, SVector{3}(CART_COORD), sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, SPH_COORD, 0.0 + 0.0im, false)
        channel_stat_signal_and_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, SPH_COORD, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, SPH_COORD, 0.0 + 0.0im, false)
        
        @test CART_COORD == sim_doa(1s, channel_stat_signal_doa_sph.enu_doa)
        @test CART_COORD == sim_doa(1s, channel_stat_interf_doa_sph.interf_enu_doa)
        @test CART_COORD == sim_doa(1s, channel_stat_signal_and_interf_doa_sph.enu_doa)
        @test CART_COORD == sim_doa(1s, channel_stat_signal_and_interf_doa_sph.interf_enu_doa)
    end

    @testset "Satellite channels with dynamic spherical DOAs" begin
        channel_dyn_signal_doa_sph = GNSSSimulator.SatelliteChannel(1, DYN_DOA_SPH, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, DYN_DOA_CART, 0.0 + 0.0im, false)
        channel_dyn_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, DYN_DOA_CART, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, DYN_DOA_SPH, 0.0 + 0.0im, false)
        channel_dyn_signal_and_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, DYN_DOA_SPH, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, DYN_DOA_SPH, 0.0 + 0.0im, false)
        
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_signal_doa_sph.enu_doa)
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_interf_doa_sph.interf_enu_doa)
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_signal_and_interf_doa_sph.enu_doa)
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_signal_and_interf_doa_sph.interf_enu_doa)
    end
end 