abstract type AbstractExistence end

struct DynamicExistence <: AbstractExistence
    existence::Vector{Bool}
    time::typeof(1.0s)
    sample_freq::typeof(1.0Hz)
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
    DynamicExistence(existence.existence, next_time, existence.sample_freq)
end

function get_existence(existence::DynamicExistence)
    get_sampled_value(existence.time, existence.sample_freq, existence.existence)
end
