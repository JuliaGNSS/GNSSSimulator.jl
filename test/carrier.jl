@testset "Carrier $type" for type in (Float32, Float64)
    carrier = StructArray{Complex{type}}(undef, 2500)
    @inferred GNSSSimulator.gen_carrier!(
        carrier,
        100.0Hz,
        2.5e6Hz,
        0.25,
        type(3.5)
    )
    @test carrier ≈ cis.(2π .* 100.0 ./ 2.5e6 .* (0:2499) .+ 2π .* 0.25) .* 3.5
end