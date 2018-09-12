@testset "Pseudo Post Correlation Signal" begin
    amplitude = 2
    amplitude_dB = 20 * log10(amplitude ^ 2) * 1dB
    ampl_std = 0.2
    phase = sqrt(pi)
    phase_std = 0.01
    static_signal = amplitude * cis(phase)
    noisy_signal_dB = GNSSSimulator.NoisyPseudoPostCorr(amplitude_dB, phase, 0.0, 0.0)
    array_noisy_complex_signals = map(signal -> sim_pseudo_post_corr_signal(3s, signal), [GNSSSimulator.NoisyPseudoPostCorr(amplitude, phase, ampl_std, phase_std) for i = 1:1000])
    noisy_signal_amplitudes = abs.(array_noisy_complex_signals)
    noisy_signal_phases = angle.(array_noisy_complex_signals)
    (mean_amplitudes, std_amplitudes) = (mean(noisy_signal_amplitudes), std(noisy_signal_amplitudes, corrected = false))
    (mean_phases, std_phases) = (mean(noisy_signal_phases), std(noisy_signal_phases, corrected = false))

    @test static_signal == sim_pseudo_post_corr_signal(3s, static_signal) 
    @test abs(static_signal) ≈ abs(sim_pseudo_post_corr_signal(3s, noisy_signal_dB))
    @test angle(static_signal) ≈ angle(sim_pseudo_post_corr_signal(3s, noisy_signal_dB))
    @test mean_amplitudes ≈ amplitude atol = 0.01
    @test std_amplitudes ≈ ampl_std atol = 0.01
    @test mean_phases ≈ phase atol = 0.001
    @test std_phases ≈ phase_std atol = 0.001
end