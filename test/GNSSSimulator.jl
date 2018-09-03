const CART_COORD = SVector{3}([1.0; 2.0; 3.0])
const SPH_COORD = SphericalFromCartesian()(CART_COORD)
const CART_COORD_VEC = [[CART_COORD]; [2 * CART_COORD]; [3 * CART_COORD]]
const DYN_DOA_CART = GNSSSimulator.DynamicDOA(CART_COORD_VEC, 1Hz)
const DYN_DOA_SPH = GNSSSimulator.DynamicDOA(SPH_COORD_VEC, 1Hz) 

@testset "Satellite channels" begin
    @testset "Channels with static spherical DOAs" begin
        channel_stat_signal_doa_sph = GNSSSimulator.SatelliteChannel(1, SPH_COORD, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, SVector{3}(CART_COORD), 0.0 + 0.0im, false)
        channel_stat_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, SVector{3}(CART_COORD), sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, SPH_COORD, 0.0 + 0.0im, false)
        channel_stat_signal_and_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, SPH_COORD, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, SPH_COORD, 0.0 + 0.0im, false)
        
        @test CART_COORD == sim_doa(1s, channel_stat_signal_doa_sph.enu_doa)
        @test CART_COORD == sim_doa(1s, channel_stat_interf_doa_sph.interf_enu_doa)
        @test CART_COORD == sim_doa(1s, channel_stat_signal_and_interf_doa_sph.enu_doa)
        @test CART_COORD == sim_doa(1s, channel_stat_signal_and_interf_doa_sph.interf_enu_doa)
    end

    @testset "Channels with dynamic spherical DOAs" begin
        channel_dyn_signal_doa_sph = GNSSSimulator.SatelliteChannel(1, DYN_DOA_SPH, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, DYN_DOA_CART, 0.0 + 0.0im, false)
        channel_dyn_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, DYN_DOA_CART, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, DYN_DOA_SPH, 0.0 + 0.0im, false)
        channel_dyn_signal_and_interf_doa_sph = GNSSSimulator.SatelliteChannel(1, DYN_DOA_SPH, sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2), true, DYN_DOA_SPH, 0.0 + 0.0im, false)
        
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_signal_doa_sph.enu_doa)
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_interf_doa_sph.interf_enu_doa)
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_signal_and_interf_doa_sph.enu_doa)
        @test CART_COORD_VEC[1] == sim_doa(0s, channel_dyn_signal_and_interf_doa_sph.interf_enu_doa)
    end
end 