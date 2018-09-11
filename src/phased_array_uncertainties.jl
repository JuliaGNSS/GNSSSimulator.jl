"""
$(SIGNATURES)

Simulates gain and phase mismatch and crosstalk. Returns a function which depends on time 't'.
"""
function sim_gain_phase_mism_and_crosstalk(
        num_ants, 
        init_crosstalk_to_direct_power,
        init_phase_mism_betw_ant_var = π / 8,
        init_gain_mism_betw_ant_var = 0.1,
        init_crosstalk_phase_var = π,
        init_crosstalk_ampl_var = init_gain_mism_betw_ant_var * uconvertp(NoUnits, init_crosstalk_to_direct_power)
    )

    init_phase_mism = randn(num_ants) * sqrt(init_phase_mism_betw_ant_var)
    init_gain_mism = ones(num_ants) + randn(num_ants) * sqrt(init_gain_mism_betw_ant_var)

    init_crosstalk_ampl = (ones(num_ants, num_ants) + randn(num_ants, num_ants) * sqrt(init_crosstalk_ampl_var)) .*
        (ones(num_ants, num_ants) - Matrix(1.0I, num_ants, num_ants)) * uconvertp(NoUnits, init_crosstalk_to_direct_power)
    init_crosstalk_phase = randn(num_ants, num_ants) * sqrt(init_crosstalk_phase_var)

    gain_and_phase_mism = init_gain_mism .* cis.(init_phase_mism)
    crosstalk = init_crosstalk_ampl .* cis.(init_crosstalk_phase)
    gain_phase_mism_and_crosstalk = normalize_gain_phase_mism_and_crosstalk(diagm(0 => gain_and_phase_mism) + crosstalk)

    t -> gain_phase_mism_and_crosstalk
end

"""
$(SIGNATURES)

Normalizes the gain and phase mismatch and crosstalk function so that the norm over columns is 1.

"""
function normalize_gain_phase_mism_and_crosstalk(gain_and_phase_mism_and_crosstalk)
    gain_and_phase_mism_and_crosstalk_norm = map(norm, julienne(gain_and_phase_mism_and_crosstalk, (:,*)))'
    gain_and_phase_mism_and_crosstalk ./ gain_and_phase_mism_and_crosstalk_norm
end