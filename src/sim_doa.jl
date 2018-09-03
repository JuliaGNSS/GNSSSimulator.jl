"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is static over time and given in Cartesian coordinates.
"""
function sim_doa(t, data::SVector)
    data
end

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is static over time and given in Spherical coordinates.
"""
function sim_doa(t, data::Spherical) 
    CartesianFromSpherical()(data)
end 

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is changing over time and given in Cartesian coordinates.
If time index exceeds data length, last available value is returned

"""
function sim_doa(t, data::DynamicDOA{Vector{T}}) where T <: SVector
    index = floor(Int, t * data.sample_freq) + 1
    index < size(data.doas, 1) ? data.doas[index] : data.doas[end]
end

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is changing over time and given in Spherical coordinates.
If time index exceeds data length, last available value is returned
"""
function sim_doa(t, data::DynamicDOA{Vector{T}}) where T <: Spherical
    index = floor(Int, t * data.sample_freq) + 1
    CartesianFromSpherical()(index < size(data.doas, 1) ? data.doas[index] : data.doas[end])
end

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference which is linearly changing over time and given in Spherical coo.
"""
function sim_doa(t, data::LinearDynamicDOA)
    CartesianFromSpherical()(Spherical(data.init_DOA.r, data.init_DOA.θ + data.Δazimuth * t, data.init_DOA.ϕ + data.Δelevation * t))
end