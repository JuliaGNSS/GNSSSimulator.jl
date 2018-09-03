#= lothars_doas = [0.6409    0.5260   -0.6634    0.8138   -0.5000   -0.9513   -0.6634         0    0.4924   -0.3100         0;
               -0.6409   -0.0646    0.3830   -0.2962   -0.5000   -0.1677   -0.5567   -0.0872    0.4132    0.8517   -0.9659;
                0.4226    0.8480    0.6428    0.5000    0.7071    0.2588    0.5000    0.9962    0.7660    0.4226    0.2588]    =#    

"""
$(SIGNATURES)

Generates an array of 4 satellite channels with default doas ("Lothar's DOAs") and no interferences.
"""
function gen_example_sat_channels()
    lothars_doas = [0.6409    0.5260   -0.6634    0.8138 ;
                    -0.6409   -0.0646    0.3830   -0.2962;
                    0.4226    0.8480    0.6428    0.5000]
    [SatelliteChannel(i, SVector{3}(lothars_doas[:,i]), NoisyPseudoPostCorr(0dB, pi/2, 0.0, 0.0), true) for i = 1:4]
end