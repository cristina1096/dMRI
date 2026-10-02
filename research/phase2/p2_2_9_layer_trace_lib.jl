"""
    sampled_layer_fraction(a, b, h; n=10_000)

Time-averaged number of layers a straight piece from x = a to x = b sits in, by brute-force sampling of
`n` midpoints (walls at multiples of W_REF, layer of thickness `h` on both faces). Deliberately independent
of MCMR's `Layers` code, so it can check it. Sampling error ≤ (number of boundary crossings)/n.
"""
function sampled_layer_fraction(a::Float64, b::Float64, h::Float64; n=10_000)
    total = 0
    for k in 1:n
        u = mod(a + (b - a) * (k - 0.5) / n, W_REF)
        total += (u <= h) + (W_REF - u <= h)
    end
    return total / n
end
