abstract type AbstractDOA end

struct DynamicDOA{T <: Vector{SVector{3, R}} where R <: Real} <: AbstractDOA
    doas::T
    time::typeof(1.0s)
    sample_freq::typeof(1.0Hz)
end

function DynamicDOA(doas::Vector{Spherical{T}}, time, sample_freq) where T
    cart_doas = sph2cart.(doas)
    DynamicDOA(cart_doas, time, sample_freq)
end

@with_kw struct LinearDynamicDOA{R <: Real} <: AbstractDOA
    doa::Spherical{R}
    Δazimuth::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
    Δelevation::typeof(1.0rad / 1.0s) = 0.0rad / 1.0s
end

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is
static over time and given in Cartesian coordinates.
"""
function propagate(doa::SVector, Δt, rng)
    doa
end

function get_doa(doa::SVector)
    doa
end


"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is
static over time and given in Spherical coordinates.
"""
function propagate(doa::Spherical, Δt, rng)
    doa
end

function get_doa(doa::Spherical)
    sph2cart(doa)
end

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference signal which is
changing over time.
If time index exceeds data length, last available value is returned
"""
function propagate(doa::DynamicDOA, Δt, rng)
    next_time = doa.time + Δt
    DynamicDOA(doa.doas, next_time, doa.sample_freq)
end

function get_doa(doa::DynamicDOA)
    get_sampled_value(doa.time, doa.sample_freq, doa.doas)
end

"""
$(SIGNATURES)

Simulates the direction of arrival of a satellite signal or interference which is linearly
changing over time and given in Spherical coordinates.
"""
function propagate(doa::LinearDynamicDOA, Δt, rng)
    next_doa = Spherical(
        doa.doa.r,
        doa.doa.θ + doa.Δazimuth * Δt,
        doa.doa.ϕ + doa.Δelevation * Δt
    )
    LinearDynamicDOA(next_doa, doa.Δazimuth, doa.Δelevation)
end
function get_doa(doa::LinearDynamicDOA)
    doa.doa
end
