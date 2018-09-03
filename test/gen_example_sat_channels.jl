@testset "Example Channels" begin
        ex_sat_channels = gen_example_sat_channels()

        @test ex_sat_channels[1].channel == 1
        @test ex_sat_channels[1].enu_doa == LOTHARS_DOAS[:,1]
        @test sim_pseudo_post_corr_signal(1s, ex_sat_channels[1].signal) == sqrt(uconvertp(NoUnits, 0dB)) * cis(pi/2)
        @test ex_sat_channels[1].exists == true
        @test ex_sat_channels[1].interf_enu_doa == SVector{3}(0.0,0.0,1.0)
        @test ex_sat_channels[1].interf_signal == 0.0 + 0.0im
        @test ex_sat_channels[1].interf_exists == false
    end