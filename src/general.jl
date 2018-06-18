"""
$(SIGNATURES)

Simulates the direction of arrivals over time `t`. This is the default configuration from Lothar Kurz.

# Examples
```julia-repl
julia> doas = sim_doas();
julia> doas(0, [trues(4); falses(7)])
3×4 Array{Float64,2}:
  0.6409   0.526   -0.6634   0.8138
 -0.6409  -0.0646   0.383   -0.2962
  0.4226   0.848    0.6428   0.5
```
"""
function sim_doas()
    doas = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
           -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
            0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]
    t -> doas
end

"""
$(SIGNATURES)

Simulates the direction over time based on the data `doas_over_time` with the sample frequency `sample_freq`.
The first dimension of the data should contain the cartesian unit vectors. The second dimension should be the time.

# Examples
```julia-repl
julia> doas_data = repeat([0.6409   0.526   -0.6634   0.8138;
                          -0.6409  -0.0646   0.383   -0.2962;
                           0.4226   0.848    0.6428   0.5   ], outer = [1,1,10]);
julia> doas = sim_doas(doas_data, 10);
julia> doas(0, trues(4))
3×4 Array{Float64,2}:
  0.6409   0.526   -0.6634   0.8138
 -0.6409  -0.0646   0.383   -0.2962
  0.4226   0.848    0.6428   0.5
```
"""
function sim_doas(doas_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        doas_over_time[:,:,index]
    end
end

"""
$(SIGNATURES)

Simulates one interference direction of arrival for all satellites.

# Examples
```julia-repl
julia> doas = sim_interf_doas(5)
julia> doas(0, trues(5))
3×5 Array{Float64,2}:
 0.664463  0.664463  0.664463  0.664463  0.664463
 0.664463  0.664463  0.664463  0.664463  0.664463
 0.34202   0.34202   0.34202   0.34202   0.34202
```
"""
function sim_interf_doas(num_sats, sph_doa::Spherical = Spherical(1.0, 45 * π / 180, 20 * π / 180))
    doas = reduce(hcat, fill(CartesianFromSpherical()(sph_doa), num_sats))
    t -> doas
end

"""
$(SIGNATURES)

Simulates the interference direction over time based on the data `doas_over_time` with the sample frequency `sample_freq`.
The first dimension of the data should contain the cartesian unit vectors. The second dimension should be the time.

# Examples
```julia-repl
julia> doas_data = repeat([0.6409   0.526   -0.6634   0.8138;
                          -0.6409  -0.0646   0.383   -0.2962;
                           0.4226   0.848    0.6428   0.5   ], outer = [1,1,10]);
julia> doas = sim_doas(doas_data, 10);
julia> doas(0, trues(4))
3×4 Array{Float64,2}:
  0.6409   0.526   -0.6634   0.8138
 -0.6409  -0.0646   0.383   -0.2962
  0.4226   0.848    0.6428   0.5
```
"""
function sim_interf_doas(doas_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        doas_over_time[:,:,index]
    end
end

"""
$(SIGNATURES)

Simulates a static satellite existence over time `t` based on the boolean array `existing_sats`. 

# Examples
```julia-repl
julia> existing_sats = sim_existing_sats(trues(4));
julia> existing_sats(0)
4-element BitArray{1}:
 true
 true
 true
 true
```
"""
function sim_existing_sats(existing_sats)
    t -> existing_sats
end

"""
$(SIGNATURES)

Simulates a varying satellite existence over time `t` based on the data `existing_sats_over_time`
and sample frequency `sample_freq`. The first dimension should hold existence of satellites and
the second dimension the time.

# Examples
```julia-repl
julia> existing_sats_data = repeat(trues(4), outer = [1,10]);
julia> existing_sats = sim_existing_sats(existing_sats_data, 10);
julia> existing_sats(0)
4-element BitArray{1}:
 true
 true
 true
 true
```
"""
function sim_existing_sats(existing_sats_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        existing_sats_over_time[:, index]
    end
end

"""
$(SIGNATURES)

Simulates a static interference existence over time `t` based on the boolean array `existing_interfs`.
Note that the length must equal to the existence of satellites. 

# Examples
```julia-repl
julia> existing_interfs = sim_existing_interfs(trues(4));
julia> existing_interfs(0)
4-element BitArray{1}:
 true
 true
 true
 true
```
"""
function sim_existing_interfs(existing_interfs)
    t -> existing_interfs
end

"""
$(SIGNATURES)

Simulates a varying interference existence over time `t` based on the data `existing_interfs_over_time`
and sample frequency `sample_freq`. The first dimension should hold existence of interferences and
the second dimension the time. Note that the length must equal to the existence of satellites.

# Examples
```julia-repl
julia> existing_sats_data = repeat(trues(4), outer = [1,10]);
julia> existing_sats = sim_existing_sats(existing_sats_data, 10);
julia> existing_sats(0)
4-element BitArray{1}:
 true
 true
 true
 true
```
"""
function sim_existing_interfs(existing_interfs_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        existing_interfs_over_time[:, index]
    end
end

"""
$(SIGNATURES)

Simulates a pseudo post correlation signal over time `t` for given `existing_sats` at that time 
instance with the power of `signal_power`. Pseudo means that no GNSS data is included.

# Examples
```julia-repl
julia> pseudo_post_corr_signal = sim_pseudo_post_corr_signal(32, 0);
julia> pseudo_post_corr_signal(0, [trues(4); falses(28)])
4-element Array{Complex{Float64},1}:
 -0.374531-0.927215im
  0.382237-0.924065im
  0.146641-0.98919im
 -0.636383+0.771374im
```
"""
function sim_pseudo_post_corr_signal(num_sats, signal_power, init_phase_var_between_signals = π)
    amplitude = sqrt(uconvertp(NoUnits, signal_power))
    init_signal_phase = randn(num_sats) * sqrt(init_phase_var_between_signals)
    signal = amplitude .* cis.(init_signal_phase)
    (t, existing_sats) -> begin
        signal .* existing_sats
    end
end

"""
$(SIGNATURES)

Simulates a varying pseudo post correlation signal over time `t` for given `existing_sats` at that time 
instance with the power of `signal_power`. Pseudo means that no GNSS data is included.

"""
function sim_pseudo_post_corr_signal(num_sats, signal_power, init_phase_var_between_signals, ampl_var, phase_var)
    pseudo_post_corr_signal = sim_pseudo_post_corr_signal(num_sats, signal_power, init_phase_var_between_signals)
    ampl_std = sqrt(ampl_var)
    phase_std = sqrt(phase_var)
    (t, existing_sats) -> begin
        signal = pseudo_post_corr_signal(t, existing_sats)
        (abs.(signal) .+ randn(num_sats) .* existing_sats .* ampl_std) .* cis.(angle.(signal) .+ randn(num_sats) .* existing_sats .* phase_std)
    end
end

"""
$(SIGNATURES)

Simulates a pseudo post correlation interference signal over time `t` for given `existing_interfs` at that time 
instance with the power of `signal_power`. Pseudo means that no GNSS data is included.

# Examples
```julia-repl
julia> pseudo_post_corr_interf_signal = sim_pseudo_post_corr_interf_signal(32, -3dB);
julia> pseudo_post_corr_interf_signal(0, [trues(4); falses(28)])
4-element Array{Complex{Float64},1}:
 -0.374531-0.927215im
  0.382237-0.924065im
  0.146641-0.98919im
 -0.636383+0.771374im
```
"""
function sim_pseudo_post_corr_interf_signal(num_interfs, signal_power, init_phase_var_between_signals = π)
    sim_pseudo_post_corr_signal(num_interfs, signal_power, init_phase_var_between_signals)
end

"""
$(SIGNATURES)

Internal pseudo post inference free signal
"""
function sim_pseudo_post_corr_interf_free_signal()
    (t, existing_sats) -> begin
        zeros(length(existing_sats))
    end
end

"""
$(SIGNATURES)

Simulates a static attitude over time `t`.

# Examples
```julia-repl
julia> attitude = sim_attitude(π / 2, π / 32, -π / 32)
julia> attitude(0)
3×3 RotXYZ{Float64}(-0.0981748, 0.0981748, 1.5708):
  6.09375e-17  -0.995185    0.0980171
  0.995185      0.00960736  0.0975452
 -0.0980171     0.0975452   0.990393
```
"""
function sim_attitude(init_yaw, init_pitch, init_roll)
    attitude = RotXYZ(init_roll, init_pitch, init_yaw)
    t -> attitude
end

"""
$(SIGNATURES)

Simulates a varying attitude over time `t` given by the data `attitude_over_time` with sample frequency
`sample_freq`. The first dimension should hold the attitude angles in the order of roll, pitch and yaw.
The second dimension should hold the time.

# Examples
```julia-repl
julia> attitude_data = repeat([-π / 32, π / 32, π / 2], outer = [1,10]);
julia> attitude = sim_attitude(attitude_data, 10);
julia> attitude(0)
3×3 RotXYZ{Float64}(-0.0981748, 0.0981748, 1.5708):
  6.09375e-17  -0.995185    0.0980171
  0.995185      0.00960736  0.0975452
 -0.0980171     0.0975452   0.990393
```
"""
function sim_attitude(attitude_over_time, sample_freq)
    t -> begin
        index = floor(Int, t * sample_freq) + 1
        RotXYZ(attitude_over_time[1, index], attitude_over_time[2, index], attitude_over_time[3, index])
    end
end

"""
$(SIGNATURES)

Simulates a varying attitude over time `t` given by the data `attitude_over_time` with sample frequency
`sample_freq` and variances given by `yaw_over_time_var`, `pitch_over_time_var` and `roll_over_time_var`. 
The first dimension should hold the attitude angles in the order of roll, pitch and yaw. The second 
dimension should hold the time.

"""
function sim_attitude(attitude_over_time, sample_freq, yaw_over_time_var, pitch_over_time_var, roll_over_time_var)
    yaw_over_time_std = sqrt(yaw_over_time_var)
    pitch_over_time_std = sqrt(pitch_over_time_var)
    roll_over_time_std = sqrt(roll_over_time_var)
    t -> begin
        yaw_noise = randn() * yaw_over_time_std
        pitch_noise = randn() * pitch_over_time_std
        roll_noise = randn() * roll_over_time_std
        index = floor(Int, t * sample_freq) + 1
        curr_att = attitude_over_time[:,index]
        RotXYZ(curr_att[1] + roll_noise, curr_att[2] + pitch_noise, curr_att[3] + yaw_noise)
    end
end

"""
$(SIGNATURES)

Simulates noise over time `t` for the `existing_sats` at that time instance with the power of
`noise_power_dB` and number of antennas `num_ants`.

# Examples
```julia-repl
julia> noise = sim_noise(-15, 4);
julia> noise(0, trues(4))
4×4 Array{Complex{Float64},2}:
  0.0758497-0.0101104im  0.0468252+0.0661756im  -0.0482784-0.229436im    0.0258948-0.116739im
 -0.0981841+0.0603515im   0.276672-0.0403435im   0.0129465-0.0112395im  -0.0486496-0.148027im
 -0.0613673-0.182071im   -0.113879+0.169017im   -0.0353773-0.111126im   -0.0873133+0.068078im
 0.00825569-0.0178933im  -0.153156+0.190351im     0.206146+0.0168782im   -0.213834-0.0693882im
```
"""
function sim_noise(noise_power, num_ants = 1)
    amplitude = sqrt(uconvertp(NoUnits, noise_power))
    (t, existing_sats) -> begin
        num_sats = length(existing_sats)
        complex.(randn(num_ants, num_sats), randn(num_ants, num_sats)) / sqrt(2) * amplitude
    end
end
