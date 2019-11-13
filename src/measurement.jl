function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters,
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {T <: Union{Float32, Float64}}
    signal = Vector{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters,
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, T <: Union{Float32, Float64}}
    signal = Matrix{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters,
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
)
    signal = Vector{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters,
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where N
    signal = Matrix{Complex{Float64}}(undef, N, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

@unroll function get_measurement!(
    signal::Union{Matrix{Complex{T}}, Vector{Complex{T}}},
    receiver::Receiver,
    emitters,
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, T <: Union{Float32, Float64}}
    if length(emitters) > 0
        phase_wraps = map(init_phase_wrap, emitters)
        C = convert_complex_or_real.(T, get_gain_phase_mism_crosstalk(receiver))
        steer_vecs = get_steer_vecs(T, emitters, receiver, manifold)
        num_samples = get_num_samples(signal)
        @inbounds @fastmath for i = 1:num_samples

            phases = calc_phases(
                emitters,
                get_intermediate_frequency(receiver),
                get_sample_frequency(receiver),
                i - 1
            )
            phase_wraps = update_phase_wraps(phase_wraps, phases, emitters)

            signal_sample = randn(rng, get_temp_signal_type(Val(N), T)) *
                get_noise_std(receiver)

            @unroll for s = 1:length(emitters)
                signal_sample = signal_sample + get_signal(
                    emitters[s],
                    phases[s],
                    phase_wraps[s],
                    steer_vecs[s],
                    rng
                )
            end

            signal_sample = C * signal_sample

            set_signal!(signal, i, signal_sample)
        end
    else
        randn!(rng, signal)
        signal .*= get_noise_std(receiver)
    end
    Δt = num_samples / get_sample_frequency(receiver)
    #next_emitters = propagate(emitters, get_intermediate_frequency(receiver), Δt, rng)
    next_receiver = propagate(receiver, Δt, rng)
    next_receiver#, next_emitters
end

function update_phase_wraps(phase_wraps, phases, emitters)
    map(phase_wraps, phases, emitters) do phase_wrap, phase, emitter
        update_phase_wrap(phase_wrap, phase, emitter)
    end
end

function calc_phases(emitters, intermediate_frequency, sample_frequency, index)
    map(emitters) do emitter
        calc_phase(emitter, index, intermediate_frequency, sample_frequency)
    end
end

function get_steer_vecs(::Type{T}, emitters, receiver, manifold) where T <: Union{Float32, Float64}
    map(emitters) do emitter
        convert_complex_or_real.(T,
            get_steer_vec(manifold, get_doa(emitter), get_attitude(receiver))
        )
    end
end

function get_temp_signal_type(::Val{N}, ::Type{T}) where {N, T}
    SVector{N, Complex{T}}
end

function get_temp_signal_type(::Val{1}, ::Type{T}) where {T}
    Complex{T}
end

function convert_complex_or_real(::Type{T}, a::Complex) where T <: Union{Float32, Float64}
    Complex{T}(a)
end

function convert_complex_or_real(::Type{T}, a) where T <: Union{Float32, Float64}
    T(a)
end

@inline Base.@propagate_inbounds function set_signal!(signals::Matrix, index, signal)
    signals[:,index] .= signal
end

@inline Base.@propagate_inbounds function set_signal!(signals::Vector, index, signal)
    signals[index] = signal
end

@inline function get_num_samples(signal::Matrix)
    size(signal, 2)
end

@inline function get_num_samples(signal::Vector)
    size(signal, 1)
end

function propagate(
    emitters::Tuple,
    intermediate_frequency,
    Δt,
    rng
)
    map(emitter -> propagate(emitter, intermediate_frequency, Δt, rng), emitters)
end
