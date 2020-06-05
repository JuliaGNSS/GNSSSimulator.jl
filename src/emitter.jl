abstract type AbstractEmitter{T} end

@inline get_doa(emitter::AbstractEmitter) = get_doa(emitter.doa)
@inline get_existence(emitter::AbstractEmitter) = get_existence(emitter.exists)

# Can be removed once https://github.com/JuliaLang/julia/pull/32968 is merged
# Julia v1.4
filteremitters(f, xs::Tuple) = Base.afoldl((ys, x) -> f(x) ? (ys..., x) : ys, (), xs...)

function get_steer_vec(manifold::AbstractManifold, emitter::AbstractEmitter, attitude)
    get_steer_vec(manifold, get_doa(emitter), attitude)
end

function get_steer_vec(manifold::IdealManifold{1}, emitter::AbstractEmitter, attitude)
    1.0
end
