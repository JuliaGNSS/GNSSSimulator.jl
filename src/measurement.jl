### Single type of emitters
function get_measurement(
    num_samples,
    receiver::AbstractReceiver{T},
    emitters::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where T <: AbstractFloat
    signal = StructArray{Complex{T}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::AbstractReceiver{T},
    emitters::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, T <: AbstractFloat}
    signal = StructArray{Complex{T}}(undef, num_samples, N)
    get_measurement!(signal, receiver, emitters, manifold, rng)
end

function get_measurement!(
    signal::Union{AbstractMatrix{Complex{T}}, AbstractVector{Complex{T}}},
    receiver::AbstractReceiver{T},
    emitters::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, T <: AbstractFloat}
    num_samples = get_num_samples(signal)
    existing_emitters = filteremitters(get_existence, emitters)
    emitter_signals = gen_signal!.(
        existing_emitters,
        get_sample_frequency(receiver),
        get_intermediate_frequency(receiver),
        get_noise_density(receiver),
        num_samples,
        Ref(rng)
    )
    steer_vecs = get_steer_vecs(T, existing_emitters, receiver, manifold)
    randn!(rng, signal)
    signal .*= get_noise_std(receiver)
    foreach(emitter_signals, steer_vecs) do emitter_signal, steer_vec
        signal .+= emitter_signal .* transpose(steer_vec)
    end
    measurement = multiply_gpmc_with_signal!(
        signal,
        get_gain_phase_mism_crosstalk(receiver)
    )

    next_receiver = propagate(receiver, num_samples, rng)
    next_emitters = propagate.(
        emitters,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        Ref(rng)
    )
    signal, next_receiver, next_emitters
end

### Two types of emitters
function get_measurement(
    num_samples,
    receiver::AbstractReceiver{T},
    emitters1::Tuple{Vararg{AbstractEmitter{T}}},
    emitters2::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where T <: AbstractFloat
    signal = StructArray{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::AbstractReceiver{T},
    emitters1::Tuple{Vararg{AbstractEmitter{T}}},
    emitters2::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, T <: AbstractFloat}
    signal = StructArray{Complex{Float64}}(undef, num_samples, N)
    get_measurement!(signal, receiver, emitters1, emitters2, manifold, rng)
end

function get_measurement!(
    signal::Union{AbstractMatrix{Complex{T}}, AbstractVector{Complex{T}}},
    receiver::AbstractReceiver{T},
    emitters1::Tuple{Vararg{AbstractEmitter{T}}},
    emitters2::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, T <: AbstractFloat}
    num_samples = get_num_samples(signal)
    existing_emitters1 = filteremitters(get_existence, emitters1)
    existing_emitters2 = filteremitters(get_existence, emitters2)
    emitter_signals1 = gen_signal!.(
        existing_emitters1,
        get_sample_frequency(receiver),
        get_intermediate_frequency(receiver),
        get_noise_density(receiver),
        num_samples,
        Ref(rng)
    )
    emitter_signals2 = gen_signal!.(
        existing_emitters2,
        get_sample_frequency(receiver),
        get_intermediate_frequency(receiver),
        get_noise_density(receiver),
        num_samples,
        Ref(rng)
    )
    steer_vecs1 = get_steer_vecs(T, existing_emitters1, receiver, manifold)
    steer_vecs2 = get_steer_vecs(T, existing_emitters2, receiver, manifold)
    randn!(rng, signal)
    signal .*= get_noise_std(receiver)
    foreach(emitter_signals1, steer_vecs1) do emitter_signal, steer_vec
        signal .+= emitter_signal .* transpose(steer_vec)
    end
    foreach(emitter_signals2, steer_vecs2) do emitter_signal, steer_vec
        signal .+= emitter_signal .* transpose(steer_vec)
    end
    measurement = multiply_gpmc_with_signal!(
        signal,
        get_gain_phase_mism_crosstalk(receiver)
    )

    next_receiver = propagate(receiver, num_samples, rng)
    next_emitters1 = propagate.(
        emitters1,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        Ref(rng)
    )
    next_emitters2 = propagate.(
        emitters2,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        Ref(rng)
    )
    measurement, next_receiver, next_emitters1, next_emitters2
end

### Three types of emitters
function get_measurement(
    num_samples,
    receiver::AbstractReceiver{T},
    emitters1::Tuple{Vararg{AbstractEmitter{T}}},
    emitters2::Tuple{Vararg{AbstractEmitter{T}}},
    emitters3::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{1} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {T <: AbstractFloat}
    signal = StructArray{Complex{Float64}}(undef, num_samples)
    get_measurement!(signal, receiver, emitters1, emitters2, emitters3, manifold, rng)
end

function get_measurement(
    num_samples,
    receiver::AbstractReceiver{T},
    emitters1::Tuple{Vararg{AbstractEmitter{T}}},
    emitters2::Tuple{Vararg{AbstractEmitter{T}}},
    emitters3::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{N},
    rng = Random.GLOBAL_RNG
) where {N, T <: AbstractFloat}
    signal = StructArray{Complex{Float64}}(undef, num_samples, N)
    get_measurement!(signal, receiver, emitters1, emitters2, emitters3, manifold, rng)
end

function get_measurement!(
    signal::Union{AbstractMatrix{Complex{T}}, AbstractVector{Complex{T}}},
    receiver::AbstractReceiver{T},
    emitters1::Tuple{Vararg{AbstractEmitter{T}}},
    emitters2::Tuple{Vararg{AbstractEmitter{T}}},
    emitters3::Tuple{Vararg{AbstractEmitter{T}}},
    manifold::AbstractManifold{N} = IdealManifold(),
    rng = Random.GLOBAL_RNG
) where {N, T <: AbstractFloat}
    num_samples = get_num_samples(signal)
    existing_emitters1 = filteremitters(get_existence, emitters1)
    existing_emitters2 = filteremitters(get_existence, emitters2)
    existing_emitters3 = filteremitters(get_existence, emitters3)
    emitter_signals1 = gen_signal!.(
        existing_emitters1,
        get_sample_frequency(receiver),
        get_intermediate_frequency(receiver),
        get_noise_density(receiver),
        num_samples,
        Ref(rng)
    )
    emitter_signals2 = gen_signal!.(
        existing_emitters2,
        get_sample_frequency(receiver),
        get_intermediate_frequency(receiver),
        get_noise_density(receiver),
        num_samples,
        Ref(rng)
    )
    emitter_signals3 = gen_signal!.(
        existing_emitters3,
        get_sample_frequency(receiver),
        get_intermediate_frequency(receiver),
        get_noise_density(receiver),
        num_samples,
        Ref(rng)
    )
    steer_vecs1 = get_steer_vecs(T, existing_emitters1, receiver, manifold)
    steer_vecs2 = get_steer_vecs(T, existing_emitters2, receiver, manifold)
    steer_vecs3 = get_steer_vecs(T, existing_emitters3, receiver, manifold)
    randn!(rng, signal)
    signal .*= get_noise_std(receiver)
    foreach(emitter_signals1, steer_vecs1) do emitter_signal, steer_vec
        signal .+= emitter_signal .* transpose(steer_vec)
    end
    foreach(emitter_signals2, steer_vecs2) do emitter_signal, steer_vec
        signal .+= emitter_signal .* transpose(steer_vec)
    end
    foreach(emitter_signals3, steer_vecs3) do emitter_signal, steer_vec
        signal .+= emitter_signal .* transpose(steer_vec)
    end
    measurement = multiply_gpmc_with_signal!(
        signal,
        get_gain_phase_mism_crosstalk(receiver)
    )

    next_receiver = propagate(receiver, num_samples, rng)
    next_emitters1 = propagate.(
        emitters1,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        Ref(rng)
    )
    next_emitters2 = propagate.(
        emitters2,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        Ref(rng)
    )
    next_emitters3 = propagate.(
        emitters3,
        num_samples,
        get_intermediate_frequency(receiver),
        get_sample_frequency(receiver),
        Ref(rng)
    )
    signal, next_receiver, next_emitters1, next_emitters2, next_emitters3
end

function convert_complex_or_real(::Type{T}, a::Complex) where T <: Union{Float32, Float64}
    Complex{T}(a)
end

function convert_complex_or_real(::Type{T}, a) where T <: Union{Float32, Float64}
    T(a)
end

function get_steer_vecs(::Type{T}, emitters, receiver, manifold) where T <: Union{Float32, Float64}
    map(emitters) do emitter
        convert_complex_or_real.(T,
            get_steer_vec(manifold, get_doa(emitter), get_attitude(receiver))
        )
    end
end

@inline function get_num_samples(signal)
    size(signal, 1)
end

function multiply_gpmc_with_signal!(signal, gpmc::Number)
    signal .*= gpmc
end

function multiply_gpmc_with_signal!(signal::AbstractMatrix{T}, gpmc::SMatrix) where T
    A_mul_Bt!(signal, gpmc)
end

function A_mul_Bt!(A::AbstractMatrix, B::SMatrix{N, N, T}) where {N, T}
    size(A, 2) == N || throw(DimensionMismatch("2nd dimension of A must have same dimension as of B"))
    @inbounds for i = 1:size(A, 1)
        v = SVector{N, T}(view(A, i, :))
        A[i, :] .= B * v
    end
    A
end

function A_mul_Bt!(A::StructArray, B::SMatrix{N, N, T}) where {N, T}
    size(A, 2) == N || throw(DimensionMismatch("2nd dimension of A must have same dimension as of B"))
    B_re = real.(B)
    B_im = imag.(B)
    @inbounds for i = 1:size(A, 1)
        v_re = SVector{N, T}(view(A.re, i, :))
        v_im = SVector{N, T}(view(A.im, i, :))
        A.re[i, :] .= B_re * v_re - B_im * v_im
        A.im[i, :] .= B_re * v_im + B_im * v_re
    end
    A
end
