
@testset "System test" begin
    sample_freq = 4e6
    interm_freq = 1_575_420_000
    gnss_system = GNSSSimulator.gpsl1_system()
    jammer_doa = GNSSSimulator.LinearDynamicDOA(init_DOA = Spherical(14, 14, 14))
    jammer = GNSSSimulator.CWJammer(1, jammer_doa, 100.0, 30.0)
    sat_doa = GNSSSimulator.LinearDynamicDOA(init_DOA = Spherical(14, 14, 14))
    sat = GNSSSimulator.Satellite(svid = 1, enu_doa = sat_doa)
    attitude = GNSSSimulator.LinearDynamicAttitude()
    sim_jammer_signal_fct = GNSSSimulator.init_sim_emitter_signal(jammer, gnss_system, sample_freq, interm_freq)
    next_sim_jammer_signal_fct, sim_jammer_signal = sim_jammer_signal_fct(4000)
    sim_sat_signal_fct = GNSSSimulator.init_sim_emitter_signal(sat, gnss_system, sample_freq, interm_freq)
    next_sim_sat_signal_fct, sim_sat_signal, sat_code_phase, sat_doppler = sim_sat_signal_fct(4000)
    emitters = [sat, jammer]
    get_steer_vec = PhasedArray.manifold(0.1904 / 4 * [1 -1 1 -1; 1 1 -1 -1; 0 0 0 0], 1575420e3)
    measurement_loop, signal = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq)

end