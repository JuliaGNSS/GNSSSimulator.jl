struct InternalStates
        doas::Matrix{Float64}
        existing_sats::Vector{Float64}
        attitude::RotXYZ
        gain_phase_mism_crosstalk::Matrix{Complex{Float64}}
        steering_vectors::Matrix{Complex{Float64}}
        signal::Vector{Complex{Float64}}
end

"""
$(SIGNATURES)

Simulates post correlation measurement. It depends on several functions:
`sim_existing_sat` which provides a boolean array of existing sats over time `t`,
`sim_post_corr_signal` which provides an array of post correlation signals over 
time `t` for current `existing_sats`, `sim_attitude` which provides the attitude
over time `t`, `sim_doas` which provides the direction of arrival over time `t` for 
current `existing_sats`, `sim_gain_phase_mism_and_crosstalk` which provides the gain
and phase mismatch and crosstalk matrix over time `t`, `sim_steering_vectors` which
provides the steering vectors over time `t` for the current `attitude` and `doas`,
and `sim_noise` which provides the noise over time `t`

"""
function sim_post_corr_measurement(
        sim_existing_sats,
        sim_post_corr_signal,
        sim_attitude,
        sim_doas,
        sim_gain_phase_mism_and_crosstalk,
        sim_steering_vectors,
        sim_noise
    )

    t -> begin
        existing_sats = sim_existing_sats(t)
        attitude = sim_attitude(t)
        doas = sim_doas(t, existing_sats)
        𝐀 = sim_steering_vectors(t, attitude, doas)
        𝐂 = sim_gain_phase_mism_and_crosstalk(t)
        𝐬 = sim_post_corr_signal(t, existing_sats)
        𝐍 = sim_noise(t, existing_sats)
        𝐘 = 𝐂 * (𝐀 .* 𝐬.' + 𝐍)
        internal_states = InternalStates(doas, existing_sats, attitude, 𝐂, 𝐀, 𝐬)
        𝐘, internal_states
    end
end

"""
$(SIGNATURES)

Simulates gain and phase mismatch and crosstalk. Returns a function which depends on time `t`.

"""
function sim_gain_phase_mism_and_crosstalk(
        num_ants, 
        init_crosstalk_to_direct_power_dB,
        init_phase_mism_betw_ant_var = π / 2,
        init_gain_mism_betw_ant_var = 0.1,
        init_crosstalk_phase_var = π,
        init_crosstalk_ampl_var = init_gain_mism_betw_ant_var * 10^(init_crosstalk_to_direct_power_dB / 10)
    )

    init_phase_mism = randn(num_ants) * sqrt(init_phase_mism_betw_ant_var)
    init_gain_mism = ones(num_ants) + randn(num_ants) * sqrt(init_gain_mism_betw_ant_var)

    init_crosstalk_ampl = (ones(num_ants, num_ants) + randn(num_ants, num_ants) * sqrt(init_crosstalk_ampl_var)) .*
        (ones(num_ants, num_ants) - eye(num_ants)) * 10^(init_crosstalk_to_direct_power_dB / 10)
    init_crosstalk_phase = randn(num_ants, num_ants) * sqrt(init_crosstalk_phase_var)

    gain_and_phase_mism = init_gain_mism .* cis.(init_phase_mism)
    crosstalk = init_crosstalk_ampl .* cis.(init_crosstalk_phase)
    gain_phase_mism_and_crosstalk = normalize_gain_phase_mism_and_crosstalk(diagm(gain_and_phase_mism) + crosstalk)

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

"""
$(SIGNATURES)

Simulates steering vectors over time using the given `get_steer_vec` function. It returns a function based on time,
attitude and doas. The output type is specified because currently the splatting make it type instable:
https://github.com/JuliaLang/julia/issues/21672

"""
function sim_steering_vectors(get_steer_vec)
    (t, attitude, doas) -> begin
        hcat(map(get_steer_vec, julienne(attitude * doas, (:,*)))...)::Array{Complex{Float64}, 2}
    end
end
