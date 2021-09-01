abstract type AbstractAmplitude end

"""
$(SIGNATURES)

Simulates the amplitude of interference signal which is
static over time.
"""
function propagate(jnr::Unitful.Gain, n0, sample_frequency, Δt, rng)
    jnr
end

@inline function get_amplitude(jnr::Unitful.Gain, n0, sample_frequency)
    sqrt(uconvertp(NoUnits, jnr) * n0 * sample_frequency)
end

@inline function get_jammer_to_noise_ratio(jnr::Unitful.Gain)
    jnr
end

struct DynamicJNR{T <: Real, VC <: AbstractVector{<: Unitful.Gain}, F <: Function} <: AbstractAmplitude
    jnrs::VC
    ampl::T
    sample_freq::typeof(1.0Hz)
    time::typeof(1.0s)
    ampl_filter::F
end

function DynamicJNR(jnrs, sample_freq; ampl_filter = x -> x, time = 0.0s)
    DynamicJNR(jnrs, 0.0, float(sample_freq), time, ampl_filter)
end

function propagate(dynamic_jnr::DynamicJNR, n0, sample_frequency, Δt, rng)
    next_time = dynamic_jnr.time + Δt
    jnr = get_sampled_value(next_time, dynamic_jnr.sample_freq, dynamic_jnr.jnrs)
    ampl = get_amplitude(jnr, n0, sample_frequency)
    ampl_filtered = dynamic_jnr.ampl_filter(ampl)
    DynamicJNR(dynamic_jnr.jnrs, ampl_filtered, dynamic_jnr.sample_freq, next_time, dynamic_jnr.ampl_filter)
end

@inline function get_amplitude(dynamic_jnr::DynamicJNR, n0, sample_frequency)
    dynamic_jnr.ampl
end

@inline function get_amplitude(dynamic_jnr::DynamicJNR)
    dynamic_jnr.ampl
end

@inline function get_jammer_to_noise_ratio(dynamic_jnr::DynamicJNR)
    dynamic_jnr
end

@inline get_amplitude(cn0::Unitful.Level, n0) = calc_amplitude_from_cn0(cn0, n0)

struct DynamicCN0{T <: Real, VC <: AbstractVector{<: Unitful.Level}, F <: Function} <: AbstractAmplitude
    cn0s::VC
    ampl::T
    sample_freq::typeof(1.0Hz)
    time::typeof(1.0s)
    ampl_filter::F
end

function DynamicCN0(cn0s, sample_freq; ampl_filter = x -> x, time = 0.0s)
    DynamicCN0(cn0s, 0.0, float(sample_freq), time, ampl_filter)
end

function propagate(dynamic_cn0::DynamicCN0, n0, Δt)
    next_time = dynamic_cn0.time + Δt
    cn0 = get_sampled_value(next_time, dynamic_cn0.sample_freq, dynamic_cn0.cn0s)
    ampl = get_amplitude(cn0, n0)
    ampl_filtered = dynamic_cn0.ampl_filter(ampl)
    DynamicCN0(dynamic_cn0.cn0s, ampl_filtered, dynamic_cn0.sample_freq, next_time, dynamic_cn0.ampl_filter)
end

@inline function get_amplitude(dynamic_cn0::DynamicCN0, n0)
    dynamic_cn0.ampl
end

@inline function get_amplitude(dynamic_cn0::DynamicCN0)
    dynamic_cn0.ampl
end

@inline function get_jammer_to_noise_ratio(dynamic_cn0::DynamicCN0)
    dynamic_cn0
end