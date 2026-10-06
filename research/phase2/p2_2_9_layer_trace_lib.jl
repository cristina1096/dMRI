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

"""
    sampled_profile_average(a, b, h; shape=:step, n=10_000)

Time-average over a straight piece x = a → b of g(d/h) summed over both faces of the gap (walls at multiples
of W_REF), by brute-force sampling of `n` midpoints. Independent of MCMR's `Layers` code. Valid for h ≤ W_REF
(only the two faces of the gap can reach it).
"""
function sampled_profile_average(a::Float64, b::Float64, h::Float64; shape=:step, n=10_000)
    g(u) = 0 <= u <= 1 ? (shape == :step ? 1.0 : 2 * (1 - u)) : 0.0
    total = 0.0
    for k in 1:n
        u = mod(a + (b - a) * (k - 0.5) / n, W_REF)
        total += g(u / h) + g((W_REF - u) / h)
    end
    return total / n
end
