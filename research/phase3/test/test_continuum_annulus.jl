@testset "annulus continuum (P3.7b Task 1)" begin
    R, D = 1.0, 3.0
    Ro25, Ro4 = krogh_radius(2.5), krogh_radius(4.0)
    @test Ro25 ≈ 2.5 / sqrt(π)

    # no sink, no bulk: S ≡ 1 up to the eigensolver roundoff floor ε‖M‖t (‖M‖ ≈ 4D/dr² ≈ 3e8 /ms at N = 2000, s = 2.5)
    let t = [0.0, 10.0, 100.0]
        S = annulus_signal(R=R, R_out=Ro25, D=D, times=t)
        @test all(abs.(S .- 1) .<= 1e-15 .+ 3e-8 .* t)
    end

    # ring integral over the full support = h R κ_out, κ_out = 1 + c̄/R; and against midpoint quadrature
    for (shape, h, L) in (("step", 0.2, 0.2), ("linear", 0.3, 0.3), ("exponential", 0.2, 0.5))
        full = outer_ring_integral(shape, h, R, L, R, R + L)
        @test full ≈ h * R * (1 + outer_centroid(shape, h, L) / R) atol = 1e-12
        n = 1_000_000
        dr = L / n
        c = L / h
        g(u) = shape == "step" ? 1.0 : shape == "linear" ? 2 * (1 - u) : exp(-u) / -expm1(-c)
        mid = sum(g(((k - 0.5) * dr) / h) * (R + (k - 0.5) * dr) for k in 1:n) * dr
        @test full ≈ mid atol = 1e-8
    end

    # uniform limit: step layer filling the whole annulus decays exactly at ρ/(hκ_out)
    let ρ = 0.02, h = Ro25 - R
        κ = 1 + outer_centroid("step", h, h) / R
        S = annulus_signal(R=R, R_out=Ro25, D=D, times=[10.0, 50.0], layer=("step", ρ, h, h))
        @test S ≈ exp.(-ρ .* [10.0, 50.0] ./ (h * κ)) atol = 3e-8 * 50     # roundoff floor ε‖M‖t
    end

    # mean rate at t = 0 equals ρ·2πR·kept / area (untruncated and truncated)
    meanrate(layer) = (Mx = annulus_operator(R=R, R_out=Ro25, D=D, layer=layer); F = eigen(Mx[1]);
                       c2 = (F.vectors' * sqrt.(Mx[2])) .^ 2 ./ sum(Mx[2]); sum(c2 .* F.values))
    area = π * (Ro25^2 - R^2)
    for (layer, truncated) in ((("step", 0.02, 0.2, 0.2), false), (("exponential", 0.02, 0.2, 0.5), true))
        kept = annulus_layer_kept(layer[1], layer[3], layer[4], R, Ro25)
        @test truncated ? kept < 1 : kept == 1.0
        @test meanrate(layer) ≈ 0.02 * 2π * R * kept / area rtol = 1e-6
    end

    # Robin: lowest rate against the exact Bessel root
    for Ro in (Ro25, Ro4), ρ in (0.007, 0.3, 30.0)
        @test annulus_rates(R=R, R_out=Ro, D=D, surface=ρ)[1] ≈ exact_rate_annulus_robin(R, Ro, D, ρ) rtol = 1e-6 atol = 1e-7   # + roundoff floor
    end

    # ceiling: Dirichlet limit, above any finite ρ
    for Ro in (Ro25, Ro4)
        cl = annulus_ceiling(R, Ro, D)
        @test cl ≈ annulus_rates(R=R, R_out=Ro, D=D, surface=1e9)[1] rtol = 1e-4
        @test cl > exact_rate_annulus_robin(R, Ro, D, 30.0)
    end

    # grid convergence, thin step layer
    let t = collect(0.0:1.0:10.0), lay = ("step", 0.3, 0.1, 0.1)
        a = annulus_signal(R=R, R_out=Ro25, D=D, times=t, layer=lay, N=2000)
        b = annulus_signal(R=R, R_out=Ro25, D=D, times=t, layer=lay, N=4000)
        @test maximum(abs.(a .- b)) <= 1e-6
    end
end
