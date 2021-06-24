function gen_carrier!(
    carrier::StructArray{<:Complex{T}},
    frequency,
    sampling_frequency,
    start_phase,
    amplitude::T
) where T
    c_re = carrier.re; c_im = carrier.im
    phase_step = upreferred(frequency / sampling_frequency)
    @avx for i in 1:length(carrier)
        c_im_temp, c_re_temp =
            sincos(T(2π) * ((i - 1) * T(phase_step) + T(start_phase)))
        c_im[i] = c_im_temp * amplitude
        c_re[i] = c_re_temp * amplitude
    end
    carrier
end