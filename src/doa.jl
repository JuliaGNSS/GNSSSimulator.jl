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
    doa::Spherical{R, R}
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

@with_kw struct RandomWalkDOA <: AbstractDOA
    azimuth::SVector{2, Float64} = SVector(0.0,0.0)
    elevation::SVector{2, Float64} = SVector(0.0,0.0)
    azimuth_vel_std::Float64
    elevation_vel_std::Float64
    elevation_limit::Float64
end

@inline function propagate(doa::RandomWalkDOA, Δt, rng)
    P = get_process(Order(2), Δt / s)
    C = cholesky(get_process_covariance(Order(2), Δt / s)).L
    elevation_limit = doa.elevation_limit
    next_azimuth = P * doa.azimuth + C * randn(rng, SVector{2,Float64}) * doa.azimuth_vel_std
    next_elevation = soft_bound(P * doa.elevation , C * randn(rng, SVector{2,Float64}) * doa.elevation_vel_std, π / 2, elevation_limit)   
    RandomWalkDOA(next_azimuth, next_elevation, doa.azimuth_vel_std, doa.elevation_vel_std, elevation_limit)
end

function get_doa(doa::RandomWalkDOA)
    Spherical(1.0, doa.azimuth[1], doa.elevation[1])
end