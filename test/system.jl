using PyPlot
#For now
function beamform(x)
    [0.5 0.5 0.5 0.5] * x
end

 @testset "System test" begin
    sample_freq = 4e6
    interm_freq = 100_000
    gnss_system = GNSSSimulator.gpsl1_system()
    jammer_doa = GNSSSimulator.LinearDynamicDOA(init_DOA = Spherical(1.0, 0.0, 0.0))
    jammer = GNSSSimulator.CWJammer(1, jammer_doa, 0.0, -50.0)
    sat_doa = GNSSSimulator.LinearDynamicDOA(init_DOA = Spherical(1.0, pi/2, pi/2))
    sat = GNSSSimulator.Satellite(svid = 1, enu_doa = sat_doa)
    attitude = GNSSSimulator.LinearDynamicAttitude()
  #   sim_jammer_signal_fct = GNSSSimulator.init_sim_emitter_signal(jammer, gnss_system, sample_freq, interm_freq)
  #  next_sim_jammer_signal_fct, sim_jammer_signal = sim_jammer_signal_fct(4000)
  #  sim_sat_signal_fct = GNSSSimulator.init_sim_emitter_signal(sat, gnss_system, sample_freq, interm_freq)
   # next_sim_sat_signal_fct, sim_sat_signal, sat_code_phase, sat_doppler = sim_sat_signal_fct(4000) 
    emitters = [sat, jammer]
    emitters = [sat]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    #get_steer_vec = sim_steering_vectors(a -> [a[1] + 0.0im, a[1] + 0.0im, a[2] + 0.0im, a[3] + 0.0im])
    measurement_loop = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq)
    next_measurement_loop, signal, internal_states = measurement_loop(4000)
    #plot(fftshift(fft(signal[1,:])))
    gen_sampled_code, get_code_phase = Tracking.init_gpsl1_codes()
    internal_code_phase = internal_states[1].code_phase
    internal_doppler = internal_states[1].doppler
    println("Internal doppler: ", internal_states[1].doppler, " internal phase: ", internal_states[1].code_phase)
    track = Tracking.init_tracking(internal_states[1].carrier_phase, 1575.42e6, interm_freq, internal_states[1].doppler, internal_states[1].code_phase, 1.023e6, sample_freq, 18.0, 1.0, 1, gen_sampled_code, get_code_phase)
    next_track, code_phase, prompts_correlated_signals, carrier_freq = track(signal, beamform)
    println("1.  ", code_phase, prompts_correlated_signals, carrier_freq-interm_freq)
    carrier_frequencies = zeros(3000)
    for i = 1:3000
    next_measurement_loop, signal, internal_states = next_measurement_loop(4000)
    next_track, code_phase, prompts_correlated_signals, carrier_freq = next_track(signal, beamform)
    carrier_frequencies[i] = carrier_freq 
    end
    figure()
    plot(carrier_frequencies)
    @test code_phase ≈ internal_states[1].code_phase atol = 3
    @test carrier_freq ≈ 50 + internal_states[1].doppler atol = 150
    println("after 1000  ", "code phase", code_phase," prompots corr sig ", prompts_correlated_signals," carrier interm freq ", carrier_freq-interm_freq, " internal code phase", internal_states[1].code_phase)
end 

@testset "continuity of measurement_loop" begin 
    sample_freq = 4e2
    interm_freq = 1_575_420_000
    gnss_system = GNSSSimulator.gpsl1_system()
    jammer_doa = GNSSSimulator.LinearDynamicDOA(init_DOA = Spherical(1.0, 0.0, 0.0))
    jammer = GNSSSimulator.CWJammer(1, jammer_doa, 0.0, -50.0)
    sat_doa = GNSSSimulator.LinearDynamicDOA(init_DOA = Spherical(1.0, pi/2, pi/2))
    sat = GNSSSimulator.Satellite(svid = 1, enu_doa = sat_doa)
    attitude = GNSSSimulator.LinearDynamicAttitude()
    emitters = [sat, jammer]
    get_steer_vec(doa) = [0.5, 0.5, 0.5, 0.5]
    measurement_loop = GNSSSimulator.init_sim_measurement(emitters, gnss_system, attitude, get_steer_vec, sample_freq, interm_freq)

    next_measurement_loop, signal_long, internal_states = measurement_loop(8)
    next_measurement_loop, signal_long2, internal_states = measurement_loop(8)
    @test signal_long == signal_long2
    next_measurement_loop, signal_short1, internal_states = measurement_loop(4)
    next_measurement_loop, signal_short2, internal_states = next_measurement_loop(4)
    short_signal = [signal_short1 signal_short2]
    @test signal_long ≈ short_signal
    @test signal_long[1,:] ≈ short_signal[1,:]
    @test signal_long[2,:] ≈ short_signal[2,:]
    @test signal_long[3,:] ≈ short_signal[3,:]
    @test signal_long[4,:] ≈ short_signal[4,:]
    #println("signal", signal_short1[40])
    #if schleife mit vergleich und i ausgeben
end
