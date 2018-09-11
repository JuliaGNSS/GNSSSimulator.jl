@testset "Attitude" begin
    @testset "Static Attitude" begin
        @test STAT_ATT == sim_attitude(3s, STAT_ATT) 

        array_struct_noisy_stat_att = [GNSSSimulator.NoisyStaticAttitude(STAT_ATT, STD_ROLL, STD_PITCH, STD_YAW) for i = 1:1000]
        meas_static_std_roll = std([sim_attitude(0s, array_struct_noisy_stat_att[i]).theta1 for i = 1:length(array_struct_noisy_stat_att)])
        meas_static_std_pitch = std([sim_attitude(0s, array_struct_noisy_stat_att[i]).theta2 for i = 1:length(array_struct_noisy_stat_att)])
        meas_static_std_yaw = std([sim_attitude(0s, array_struct_noisy_stat_att[i]).theta3 for i = 1:length(array_struct_noisy_stat_att)])

        @test meas_static_std_roll ≈ STD_ROLL atol = 0.001
        @test meas_static_std_pitch ≈ STD_PITCH atol = 0.001
        @test meas_static_std_yaw ≈ STD_YAW atol = 0.003
    end

    @testset "Dynamic Attitude" begin
        struct_dyn_att = GNSSSimulator.DynamicAttitude(DYN_ATT, 1Hz)
        struct_lin_dyn_att = GNSSSimulator.LinearDynamicAttitude(STAT_ATT, ΔROLL_PER_S, ΔPITCH_PER_S, ΔYAW_PER_S)

        @test DYN_ATT[3] == sim_attitude(2s, struct_dyn_att)
        @test DYN_ATT[3] == sim_attitude(5s, struct_dyn_att) # time exceeds data length --> return last available value 
        @test DYN_ATT[3] ≈ sim_attitude(2s, struct_lin_dyn_att)

        array_struct_noisy_dyn_att = [GNSSSimulator.NoisyDynamicAttitude(DYN_ATT, 1Hz, STD_ROLL, STD_PITCH, STD_YAW) for i = 1:1000]
        meas_dyn_std_roll = std([sim_attitude(2s, array_struct_noisy_dyn_att[i]).theta1 for i = 1:length(array_struct_noisy_dyn_att)])
        meas_dyn_std_pitch = std([sim_attitude(2s, array_struct_noisy_dyn_att[i]).theta2 for i = 1:length(array_struct_noisy_dyn_att)])
        meas_dyn_std_yaw = std([sim_attitude(2s, array_struct_noisy_dyn_att[i]).theta3 for i = 1:length(array_struct_noisy_dyn_att)])

        array_struct_lin_noisy_dyn_att = [GNSSSimulator.NoisyLinearDynamicAttitude(STAT_ATT, ΔROLL_PER_S, ΔPITCH_PER_S, ΔYAW_PER_S, STD_ROLL, STD_PITCH, STD_YAW) for i = 1:1000]
        meas_lin_dyn_std_roll = std([sim_attitude(2s, array_struct_noisy_dyn_att[i]).theta1 for i = 1:length(array_struct_noisy_dyn_att)])
        meas_lin_dyn_std_pitch = std([sim_attitude(2s, array_struct_noisy_dyn_att[i]).theta2 for i = 1:length(array_struct_noisy_dyn_att)])
        meas_lin_dyn_std_yaw = std([sim_attitude(2s, array_struct_noisy_dyn_att[i]).theta3 for i = 1:length(array_struct_noisy_dyn_att)]) 

        @test meas_dyn_std_roll ≈ STD_ROLL atol = 0.002
        @test meas_dyn_std_pitch ≈ STD_PITCH atol = 0.002
        @test meas_dyn_std_yaw ≈ STD_YAW atol = 0.002
        @test meas_lin_dyn_std_roll ≈ STD_ROLL atol = 0.002
        @test meas_lin_dyn_std_pitch ≈ STD_PITCH atol = 0.002
        @test meas_lin_dyn_std_yaw ≈ STD_YAW atol = 0.002
    end
end