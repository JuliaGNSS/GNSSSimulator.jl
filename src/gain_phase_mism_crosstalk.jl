abstract type AbstractGainPhaseMismCrosstalk end

function propagate(gain_phase_mism_crosstalk::SMatrix, Δt)
    gain_phase_mism_crosstalk
end

function get_gain_phase_mism_crosstalk(gain_phase_mism_crosstalk)
    gain_phase_mism_crosstalk
end

"""
$(SIGNATURES)

Normalizes the gain and phase mismatch and crosstalk function so that the norm over columns is 1.

"""
function normalize_gain_phase_mism_and_crosstalk(gain_and_phase_mism_and_crosstalk)
    gain_and_phase_mism_and_crosstalk_norm = map(norm, Slices(gain_and_phase_mism_and_crosstalk, False(), True()))'
    gain_and_phase_mism_and_crosstalk ./ gain_and_phase_mism_and_crosstalk_norm
end
