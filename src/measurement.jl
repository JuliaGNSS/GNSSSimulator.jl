### Single type of emitters
function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters::NTuple{NS, <:AbstractEmitter},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {NS, T <: Union{Float32, Float64}}
    signal = Vector{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters::NTuple{NS, <:AbstractEmitter},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, NS, T <: Union{Float32, Float64}}
    signal = Matrix{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters::NTuple{NS, <:AbstractEmitter},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {NS}
    signal = Vector{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters::NTuple{NS, <:AbstractEmitter},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, NS}
    signal = Matrix{Complex{Float64}}(undef, N, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

@unroll function get_measurement!(
    signal::Union{AbstractMatrix{Complex{T}}, AbstractVector{Complex{T}}},
    receiver::Receiver,
    emitters::NTuple{NS, <:AbstractEmitter},
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, NS, T <: Union{Float32, Float64}}
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
    next_receiver = propagate(receiver, num_samples, rng)
    next_emitters = propagate(
        emitters,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        rng
    )
    signal, next_receiver, next_emitters
end

### Two types of emitters
function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {NS1, NS2, T <: Union{Float32, Float64}}
    signal = Vector{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, manifold, rng)
end

function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, NS1, NS2, T <: Union{Float32, Float64}}
    signal = Matrix{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {NS1, NS2}
    signal = Vector{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, NS1, NS2}
    signal = Matrix{Complex{Float64}}(undef, N, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, manifold, rng)
end

@unroll function get_measurement!(
    signal::Union{AbstractMatrix{Complex{T}}, AbstractVector{Complex{T}}},
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, NS1, NS2, T <: Union{Float32, Float64}}
    phase_wraps1 = map(init_phase_wrap, emitters1)
    phase_wraps2 = map(init_phase_wrap, emitters2)
    C = convert_complex_or_real.(T, get_gain_phase_mism_crosstalk(receiver))
    steer_vecs1 = get_steer_vecs(T, emitters1, receiver, manifold)
    steer_vecs2 = get_steer_vecs(T, emitters2, receiver, manifold)
    num_samples = get_num_samples(signal)
    @inbounds @fastmath for i = 1:num_samples

        phases1 = calc_phases(
            emitters1,
            get_intermediate_frequency(receiver),
            get_sample_frequency(receiver),
            i - 1
        )
        phase_wraps1 = update_phase_wraps(phase_wraps1, phases1, emitters1)
        phases2 = calc_phases(
            emitters2,
            get_intermediate_frequency(receiver),
            get_sample_frequency(receiver),
            i - 1
        )
        phase_wraps2 = update_phase_wraps(phase_wraps2, phases2, emitters2)

        signal_sample = randn(rng, get_temp_signal_type(Val(N), T)) *
            get_noise_std(receiver)

        @unroll for s = 1:length(emitters1)
            signal_sample = signal_sample + get_signal(
                emitters1[s],
                phases1[s],
                phase_wraps1[s],
                steer_vecs1[s],
                rng
            )
        end
        @unroll for s = 1:length(emitters2)
            signal_sample = signal_sample + get_signal(
                emitters2[s],
                phases2[s],
                phase_wraps2[s],
                steer_vecs2[s],
                rng
            )
        end

        signal_sample = C * signal_sample
        set_signal!(signal, i, signal_sample)
    end
    next_receiver = propagate(receiver, num_samples, rng)
    next_emitters1 = propagate(
        emitters1,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        rng
    )
    next_emitters2 = propagate(
        emitters2,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        rng
    )
    signal, next_receiver, next_emitters1, next_emitters2
end

### Three types of emitters
function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    emitters3::NTuple{NS3, <:AbstractEmitter},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {NS1, NS2, NS3, T <: Union{Float32, Float64}}
    signal = Vector{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, emitters3, manifold, rng)
end

function get_measurement(
    ::Type{T},
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    emitters3::NTuple{NS3, <:AbstractEmitter},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, NS1, NS2, NS3, T <: Union{Float32, Float64}}
    signal = Matrix{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, emitters3, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    emitters3::NTuple{NS3, <:AbstractEmitter},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {NS1, NS2, NS3}
    signal = Vector{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, emitters3, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    emitters3::NTuple{NS3, <:AbstractEmitter},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, NS1, NS2, NS3}
    signal = Matrix{Complex{Float64}}(undef, N, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, emitters3, manifold, rng)
end

@unroll function get_measurement!(
    signal::Union{AbstractMatrix{Complex{T}}, AbstractVector{Complex{T}}},
    receiver::Receiver,
    emitters1::NTuple{NS1, <:AbstractEmitter},
    emitters2::NTuple{NS2, <:AbstractEmitter},
    emitters3::NTuple{NS3, <:AbstractEmitter},
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, NS1, NS2, NS3, T <: Union{Float32, Float64}}
    phase_wraps1 = map(init_phase_wrap, emitters1)
    phase_wraps2 = map(init_phase_wrap, emitters2)
    phase_wraps3 = map(init_phase_wrap, emitters3)
    C = convert_complex_or_real.(T, get_gain_phase_mism_crosstalk(receiver))
    steer_vecs1 = get_steer_vecs(T, emitters1, receiver, manifold)
    steer_vecs2 = get_steer_vecs(T, emitters2, receiver, manifold)
    steer_vecs3 = get_steer_vecs(T, emitters3, receiver, manifold)
    num_samples = get_num_samples(signal)
    @inbounds @fastmath for i = 1:num_samples

        phases1 = calc_phases(
            emitters1,
            get_intermediate_frequency(receiver),
            get_sample_frequency(receiver),
            i - 1
        )
        phase_wraps1 = update_phase_wraps(phase_wraps1, phases1, emitters1)
        phases2 = calc_phases(
            emitters2,
            get_intermediate_frequency(receiver),
            get_sample_frequency(receiver),
            i - 1
        )
        phase_wraps2 = update_phase_wraps(phase_wraps2, phases2, emitters2)
        phases3 = calc_phases(
            emitters3,
            get_intermediate_frequency(receiver),
            get_sample_frequency(receiver),
            i - 1
        )
        phase_wraps3 = update_phase_wraps(phase_wraps3, phases3, emitters3)

        signal_sample = randn(rng, get_temp_signal_type(Val(N), T)) *
            get_noise_std(receiver)

        @unroll for s = 1:length(emitters1)
            signal_sample = signal_sample + get_signal(
                emitters1[s],
                phases1[s],
                phase_wraps1[s],
                steer_vecs1[s],
                rng
            )
        end
        @unroll for s = 1:length(emitters2)
            signal_sample = signal_sample + get_signal(
                emitters2[s],
                phases2[s],
                phase_wraps2[s],
                steer_vecs2[s],
                rng
            )
        end
        @unroll for s = 1:length(emitters3)
            signal_sample = signal_sample + get_signal(
                emitters3[s],
                phases3[s],
                phase_wraps3[s],
                steer_vecs3[s],
                rng
            )
        end

        signal_sample = C * signal_sample
        set_signal!(signal, i, signal_sample)
    end
    next_receiver = propagate(receiver, num_samples, rng)
    next_emitters1 = propagate(
        emitters1,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        rng
    )
    next_emitters2 = propagate(
        emitters2,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        rng
    )
    next_emitters3 = propagate(
        emitters3,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        rng
    )
    signal, next_receiver, next_emitters1, next_emitters2, next_emitters3
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

@inline Base.@propagate_inbounds function set_signal!(
    signals::AbstractMatrix,
    index,
    signal
)
    signals[:,index] .= signal
end

@inline Base.@propagate_inbounds function set_signal!(
    signals::AbstractVector,
    index,
    signal
)
    signals[index] = signal
end

@inline function get_num_samples(signal::AbstractMatrix)
    size(signal, 2)
end

@inline function get_num_samples(signal::AbstractVector)
    size(signal, 1)
end

function propagate(
    emitters::NTuple{N, <: AbstractEmitter},
    num_samples::Integer,
    intermediate_frequency,
    sample_frequency,
    rng
) where N
    map(emitters) do emitter
        propagate(emitter, num_samples, intermediate_frequency, sample_frequency, rng)
    end
end
