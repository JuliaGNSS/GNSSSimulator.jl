abstract type AbstractEmitter end
abstract type AbstractEmitterPhases end

@inline get_amplitude(emitter::AbstractEmitter) = emitter.amplitude
@inline get_doa(emitter::AbstractEmitter) = get_doa(emitter.doa)
@inline get_existence(emitter::AbstractEmitter) = get_existence(emitter.exists)
