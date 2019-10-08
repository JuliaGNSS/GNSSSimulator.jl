abstract type AbstractExistence end

struct DynamicExistence <: AbstractExistence
    existence::Vector{Bool}
    sample_freq::typeof(1.0Hz)
    time::typeof(1.0s)
end

function DynamicExistence(existence, sample_freq; time = 0.0s)
    DynamicExistence(existence, sample_freq, time)
end

"""
$(SIGNATURES)

Simulates the static existence of either a satellite signal or a interference signal. Type: Boolean (true if signal exists).
"""
function propagate(existence::Bool, Δt, rng)
    existence
end

function get_existence(existence::Bool)
    existence
end

"""
$(SIGNATURES)

Simulates the dynamic existence of either a satellite signal or a interference signal. Type: Boolean (true if signal exists).
If time index exceeds data length, last available value is returned
"""
function propagate(existence::DynamicExistence, Δt, rng)
    next_time = existence.time + Δt
    DynamicExistence(existence.existence, existence.sample_freq, next_time)
end

function get_existence(existence::DynamicExistence)
    get_sampled_value(existence.time, existence.sample_freq, existence.existence)
end
