"""
$(SIGNATURES)

Simulates the static existence of either a satellite signal or a interference signal. Type: Boolean ('1' if signal exists).
"""
function sim_existence(t, existence::Bool)
    existence
end

"""
$(SIGNATURES)

Simulates the dynamic existence of either a satellite signal or a interference signal. Type: Boolean ('1' if signal exists).
If time index exceeds data length, last available value is returned
"""
function sim_existence(t, data::DynamicExistence)
    index = floor(Int, t * data.sample_freq) + 1
    index < size(data.existence, 1) ? data.existence[index] : data.existence[end]
end