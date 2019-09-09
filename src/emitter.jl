abstract type AbstractEmitter end

get_amplitude(emitter::AbstractEmitter) = emitter.amplitude
get_doa(emitter::AbstractEmitter) = get_doa(emitter.doa)
get_existence(emitter::AbstractEmitter) = get_existence(emitter.exists)
