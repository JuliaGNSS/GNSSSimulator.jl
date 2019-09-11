struct Noise{N}
    noise::N
    noise_std::Float64
end

function propagate(noise::Noise{T}, Δt) where T
    Noise(randn(T) * noise.noise_std, noise.noise_std)
end

function get_noise(noise)
    noise.noise
end
