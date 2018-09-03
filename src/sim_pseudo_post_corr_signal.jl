"""
$(SIGNATURES)

Creates a pseudo post correlation signal for either a satellite signal or an interference signal. "Pseudo" means that no GNSS data is included.
"""
function sim_pseudo_post_corr_signal(t, signal::Complex{Float64})
    signal
end

"""
$(SIGNATURES)

Creates a pseudo post correlation signal for either a noisy satellite signal or a noisy interference signal. "Pseudo" means that no GNSS data is included.
"""
function sim_pseudo_post_corr_signal(t, signal::NoisyPseudoPostCorr)
    (signal.ampl + signal.ampl_std * randn()) .* cis(signal.phase + signal.phase_std * randn())
end