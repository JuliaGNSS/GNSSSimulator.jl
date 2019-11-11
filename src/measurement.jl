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
) where N
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

function get_measurement!(
    signal::Union{Matrix{Complex{T}}, Vector{Complex{T}}},
    receiver::Receiver,
    emitters,
    manifold = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {T <: Union{Float32, Float64}}
    if length(emitters) > 0
        phase_wraps = map(init_phase_wrap, emitters)
        steer_vecs = map(emitters) do emitter
            convert_complex_or_real.(
                T,
                get_steer_vec(manifold, get_doa(emitter), get_attitude(receiver))
            )
        end
        C = convert_complex_or_real.(T, get_gain_phase_mism_crosstalk(receiver))

        @inbounds @fastmath for i = 1:get_num_samples(signal)

            phases = map(emitters) do emitter
                calc_phase(
                    emitter,
                    i - 1,
                    get_intermediate_frequency(receiver),
                    get_sample_frequency(receiver)
                )
            end
            phase_wraps = map(phase_wraps, phases, emitters) do phase_wrap, phase, emitter
                update_phase_wrap(phase_wrap, phase, emitter)
            end
#            temp = zero(typeof(steer_vecs[1]))
#            temp = get_signal(phases[1], emitters[1], steer_vecs[1], rng)
#            for s = 1:length(emitters)
#                temp = temp + get_signal(phases[s], emitters[s], steer_vecs[s], rng)
#            end
            emitter_signals = map(
                phases,
                phase_wraps,
                emitters,
                steer_vecs
            ) do phase, phase_wrap, emitter, steer_vec
                get_signal(phase, phase_wrap, emitter, steer_vec, rng)
            end
            emitter_signal = sum(emitter_signals)
#            temp = Base.mapfoldl(pes -> get_signal(pes..., rng), +, zip(phases, emitters, steer_vecs))
#            temp = sum_over_signals(phases, emitters, steer_vecs, rng)
            #temp = sum(temps)
            emitter_signal = C * (emitter_signal + randn(rng, typeof(emitter_signal)) *
                get_noise_std(receiver))
            set_signal!(signal, i, emitter_signal)
        end
    else
        randn!(rng, signal)
        signal .*= get_noise_std(receiver)
    end
    signal
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
    emitters::Tuple{<:AbstractEmitter},
    intermediate_frequency,
    Δt,
    rng = Random.GLOBAL_RNG
)
    map(emitter -> propagate(emitter, intermediate_frequency, Δt, rng), emitters)
end
