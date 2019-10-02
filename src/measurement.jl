struct ReceivedSignal{
    ES <: Tuple,
    R <: AbstractReceiver
}
    emitters::ES
    receiver::R
end


function get_measurement(::Type{T}, num_samples, received_signal::ReceivedSignal, manifold::AbstractManifold{1} = IdealManifold(), rng = Random.GLOBAL_RNG) where {T <: Union{Float32, Float64}}
    signal = Vector{Complex{T}}(undef, num_samples)
    get_measurement!(signal, received_signal, manifold, rng)
end

function get_measurement(::Type{T}, num_samples, received_signal::ReceivedSignal, manifold::AbstractManifold{N}, rng = Random.GLOBAL_RNG) where {N, T <: Union{Float32, Float64}}
    signal = Matrix{Complex{T}}(undef, num_samples)
    get_measurement!(signal, received_signal, manifold, rng)
end

function get_measurement(num_samples, received_signal::ReceivedSignal, manifold::AbstractManifold{1} = IdealManifold(), rng = Random.GLOBAL_RNG) where N
    signal = Vector{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, received_signal, manifold, rng)
end

function get_measurement(num_samples, received_signal::ReceivedSignal, manifold::AbstractManifold{N}, rng = Random.GLOBAL_RNG) where N
    signal = Matrix{Complex{Float64}}(undef, N, num_samples)
    get_measurement!(signal, received_signal, manifold, rng)
end

function get_measurement!(signal::Union{Matrix{Complex{T}}, Vector{Complex{T}}}, received_signal::ReceivedSignal, manifold = IdealManifold(), rng = Random.GLOBAL_RNG) where {T <: Union{Float32, Float64}}
    emitters = get_emitters(received_signal)
    if length(emitters) > 0
        receiver = get_receiver(received_signal)
        Δt_fast = 1 / get_sample_frequency(receiver)
        phases = get_phase.(emitters)
        steer_vecs = map(emitter -> convert_complex_or_real.(T, get_steer_vec(manifold, get_doa(emitter), get_attitude(receiver))), emitters)
        C = convert_complex_or_real.(T, get_gain_phase_mism_crosstalk(receiver))
        @inbounds @fastmath for i = 1:get_num_samples(signal)
            temp = get_signal(phases[1], emitters[1], steer_vecs[1], rng)
            for s = 2:length(emitters)
                temp = temp + get_signal(phases[s], emitters[s], steer_vecs[s], rng)
            end
            temp = C * (temp + randn(rng, typeof(temp)) * get_noise_std(receiver))
            set_signal!(signal, i, temp)
            phases = fast_propagate.(phases, emitters, get_intermediate_frequency(receiver), Δt_fast)
        end
    else
        randn!(rng, signal)
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

@inline function get_receiver(received_signal::ReceivedSignal)
    received_signal.receiver
end

@inline function get_emitters(received_signal::ReceivedSignal)
    received_signal.emitters
end

function propagate(received_signal::ReceivedSignal, Δt, rng = Random.GLOBAL_RNG)
    emitters = map(received_signal.emitters) do emitter
        propagate(emitter, get_intermediate_frequency(get_receiver(received_signal)), Δt, rng)
    end
    receiver = propagate(received_signal.receiver, Δt, rng)
    ReceivedSignal(emitters, receiver)
end
