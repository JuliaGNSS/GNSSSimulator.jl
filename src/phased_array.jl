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
`existing_sat` which provides a boolean array of existing sats over time `t`,
`post_corr_signal` which provides an array of post correlation signals over 
time `t` for current `existing_sats`, `attitude` which provides the attitude
over time `t`, `doas` which provides the direction of arrival over time `t` for 
current `existing_sats`, `gain_phase_mism_and_crosstalk` which provides the gain
and phase mismatch and crosstalk matrix over time `t`, `steering_vectors` which
provides the steering vectors over time `t` for the current `attitude` and `doas`,
and `noise` which provides the noise over time `t`

# Examples
```julia-repl
julia> doas = sim_doas()
julia> existing_sats = sim_existing_sats(trues(11))
julia> pseudo_post_corr_signal = sim_pseudo_post_corr_signal(11, 0dB)
julia> attitude = sim_attitude(0.0, 0.0, 0.0)
julia> gain_phase_mism_and_crosstalk = sim_gain_phase_mism_and_crosstalk(4, -15dB)
julia> steering_vectors = sim_steering_vectors(a -> [a[1] + 0.0im, a[1] + 0.0im, a[2] + 0.0im, a[3] + 0.0im])
julia> noise = sim_noise(-15dB, 4)
julia> measurement = sim_post_corr_measurement(
        existing_sats,
        pseudo_post_corr_signal,
        attitude,
        doas,
        gain_phase_mism_and_crosstalk,
        steering_vectors,
        noise)
```
"""
function sim_post_corr_measurement(
        existing_sats,
        post_corr_signal,
        attitude,
        doas,
        gain_phase_mism_and_crosstalk,
        steering_vectors,
        noise
    )

    t -> begin
        curr_existing_sats = existing_sats(t)
        curr_attitude = attitude(t)
        curr_doas = doas(t, curr_existing_sats)
        𝐀 = steering_vectors(t, curr_attitude, curr_doas)
        𝐂 = gain_phase_mism_and_crosstalk(t)
        𝐬 = post_corr_signal(t, curr_existing_sats)
        𝐍 = noise(t, curr_existing_sats)
        𝐘 = 𝐂 * (𝐀 .* 𝐬.' + 𝐍)
        internal_states = InternalStates(curr_doas, curr_existing_sats, curr_attitude, 𝐂, 𝐀, 𝐬)
        𝐘, internal_states
    end
end

"""
$(SIGNATURES)

Simulates gain and phase mismatch and crosstalk. Returns a function which depends on time `t`.

# Examples
```julia-repl
julia> gain_phase_mism_and_crosstalk = sim_gain_phase_mism_and_crosstalk(4, -15);
julia> gain_phase_mism_and_crosstalk(0)
4×4 Array{Complex{Float64},2}:
  0.643835+0.764009im    0.0260979-0.000931496im   0.021727-0.0304959im    -0.011296-0.0253616im
  0.012092+0.0213692im    0.331865+0.942271im     0.0108977+0.0328915im    0.0130248+0.0264729im
   0.01644-0.0194899im  0.00935026+0.0230717im     0.987765+0.143705im   -8.69255e-5+0.0299151im
 0.0206548-0.0093261im  -0.0262344+0.00138109im   0.0249238+0.0211023im     0.883451+0.465808im
```
"""
function sim_gain_phase_mism_and_crosstalk(
        num_ants, 
        init_crosstalk_to_direct_power,
        init_phase_mism_betw_ant_var = π / 2,
        init_gain_mism_betw_ant_var = 0.1,
        init_crosstalk_phase_var = π,
        init_crosstalk_ampl_var = init_gain_mism_betw_ant_var * uconvertp(NoUnits, init_crosstalk_to_direct_power)
    )

    init_phase_mism = randn(num_ants) * sqrt(init_phase_mism_betw_ant_var)
    init_gain_mism = ones(num_ants) + randn(num_ants) * sqrt(init_gain_mism_betw_ant_var)

    init_crosstalk_ampl = (ones(num_ants, num_ants) + randn(num_ants, num_ants) * sqrt(init_crosstalk_ampl_var)) .*
        (ones(num_ants, num_ants) - eye(num_ants)) * uconvertp(NoUnits, init_crosstalk_to_direct_power)
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

# Examples
```julia-repl
julia> steering_vectors = sim_steering_vectors(a -> complex.(randn(4), randn(4)));
julia> steering_vectors(0, RotXYZ(0,0,0), [zeros(2,5); ones(1,5)])
4×5 Array{Complex{Float64},2}:
  1.21626+0.501292im  -0.081102-0.278121im     -0.8059+1.00751im   0.559731+0.638522im  0.598888+1.12544im
  1.43552-1.20201im     1.19489+0.0726409im  -0.324954-0.652084im  0.282833-0.28839im   -1.80852-1.30904im
 0.647119+0.717731im   0.573429-0.298508im    -1.07383+1.89825im    1.06957-0.317897im  0.160637+0.0455821im
  1.64984+0.230963im  0.0961698-0.45558im      1.75346-3.2289im    -0.63785+0.460729im  0.567883-0.300519im
```
"""
function sim_steering_vectors(get_steer_vec)
    (t, attitude, doas) -> begin
        reduce(hcat, map(get_steer_vec, julienne(attitude * doas, (:,*))))::Array{Complex{Float64}, 2}
    end
end
