"""
$(SIGNATURES)

Generates an array of 4 satellite channels with default doas ("Lothar's DOAs") and no interferences.
"""
function gen_example_sat_channels()
    [SatelliteChannel(i, SVector{3}(LOTHARS_DOAS[:,i]), NoisyPseudoPostCorr(0dB, pi/2, 0.0, 0.0), true) for i = 1:4]
end