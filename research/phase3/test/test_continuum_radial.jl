@testset "radial continuum: exact cases" begin
    R, D, ρ = 1.0, 3.0, 0.02
    ts = [0.0, 10.0, 50.0]
    # eigenvalue roundoff ≈ ε·‖M‖ ≈ ε·4D N²/R² ≈ 1e-8 /ms (N = 2000), so exact cases agree to ≈ t·1e-8 ≤ 5e-7 at t = 50
    tol = 2e-6

    @testset "free diffusion and bulk only" begin
        @test all(isapprox.(radial_signal(R=R, D=D, times=ts), 1.0; atol=tol))
        @test all(isapprox.(radial_signal(R=R, D=D, times=ts, R2_bulk=1 / 80), exp.(-ts ./ 80); rtol=tol))
    end

    @testset "κ by quadrature (D3.1)" begin
        @test radial_kappa("step", 0.2, 1.0) ≈ 0.9 rtol=1e-8
        @test radial_kappa("step", 1.0, 1.0) ≈ 0.5 rtol=1e-8
        @test radial_kappa("linear", 0.3, 1.0) ≈ 1 - 0.1 rtol=1e-8      # c̄ = h/3
    end

    @testset "Review Focus 2: h = R gives the exact uniform rate 2ρ/R" begin
        @test all(isapprox.(radial_signal(R=R, D=D, times=ts, layer=("step", ρ, R)), exp.(-2ρ / R .* ts); rtol=tol))
    end

    @testset "initial slope −ρ·(S/V) = −2ρ/R for every model (uniform start)" begin
        δ = 1e-6
        for kw in ((layer=("step", ρ, 0.2),), (layer=("linear", ρ, 0.3),), (surface=ρ,))
            S = radial_signal(; R=R, D=D, times=[δ], kw...)[1]
            @test (S - 1) / δ ≈ -2ρ / R rtol=1e-3
        end
    end

    @testset "lowest rate: conventional (Robin) vs D k² with D k J₁(kR) = ρ J₀(kR)" begin
        for (Rr, Dd, ρρ) in ((1.0, 3.0, 0.02), (0.5, 1.0, 0.1), (2.0, 3.0, 0.5))
            @test radial_rates(R=Rr, D=Dd, surface=ρρ, N=4000)[1] ≈ exact_rate_conventional(Rr, Dd, ρρ) rtol=1e-3
        end
    end

    @testset "lowest rate: step layer vs exact Bessel J/I/K matching (Review Focus 3: strong ρ)" begin
        for (Rr, Dd, ρρ, h) in ((1.0, 3.0, 0.02, 0.2), (1.0, 3.0, 0.5, 0.1), (0.5, 1.0, 0.05, 0.25))
            @test radial_rates(R=Rr, D=Dd, layer=("step", ρρ, h), N=2000)[1] ≈ exact_rate_step(Rr, Dd, ρρ, h) rtol=1e-5
        end
    end

    @testset "Review Focus 1: Jensen bounds, including h below the grid spacing" begin
        E = 2ρ / R * 50
        for h in (1e-4, 0.01, 0.2, 0.7)
            S = radial_signal(R=R, D=D, times=[50.0], layer=("step", ρ, h), N=1000)[1]     # dr = 1e-3
            @test isfinite(S)
            @test S >= exp(-E) * (1 - 1e-9)                       # J ≥ 0 (Jensen, a theorem; the ordering vs the conventional model is a result, not tested)
        end
    end

    @testset "grid convergence" begin
        a = radial_signal(R=R, D=D, times=[50.0], layer=("step", ρ, 0.2), N=1000)[1]
        b = radial_signal(R=R, D=D, times=[50.0], layer=("step", ρ, 0.2), N=4000)[1]
        @test abs(a / b - 1) < 1e-6
        a = radial_signal(R=R, D=D, times=[50.0], surface=ρ, N=1000)[1]
        b = radial_signal(R=R, D=D, times=[50.0], surface=ρ, N=4000)[1]
        @test abs(a / b - 1) < 1e-4
    end
end

@testset "radial continuum: Robin condition on the face is second order (final review I3)" begin
    R, D, ρ = 1.0, 3.0, 0.02
    a = radial_signal(R=R, D=D, times=[50.0], surface=ρ, N=1000)[1]
    b = radial_signal(R=R, D=D, times=[50.0], surface=ρ, N=2000)[1]
    @test abs(a - b) < 2e-8                      # cell-centred sink: 4.5e-7
    @test radial_rates(R=R, D=D, surface=ρ, N=1000)[1] ≈ exact_rate_conventional(R, D, ρ) rtol=1e-6
end

@testset "radial continuum: linear and exponential (exact ring averages)" begin
    R, D, ρ = 1.0, 3.0, 0.02
    @testset "∫ΔR2 dA = 2πRρ for every profile, including h below the grid spacing" begin
        for (shape, h) in (("linear", 1e-4), ("linear", 0.3), ("exponential", 1e-4), ("exponential", 0.2), ("exponential", 3.0))
            (M, A) = radial_operator(R=R, D=D, layer=(shape, ρ, h), N=1000)
            (M0, _) = radial_operator(R=R, D=D, N=1000)
            V = (M.dv .- M0.dv)            # diagonal difference = V (the stiffness part is identical)
            @test sum(A .* V) ≈ 2π * R * ρ rtol=1e-8      # V from a difference of diagonals ≈ D/dr² ≈ 3e6: cancellation ≈ 1e-10 per ring
        end
    end
    @testset "κ for the exponential with support R" begin
        c = R / 0.2
        @test radial_kappa("exponential", 0.2, R) ≈ 1 - 0.2 * (1 - (c + 1) * exp(-c)) / (1 - exp(-c)) / R rtol=1e-8
    end
    @testset "initial slope −2ρ/R and Jensen bound" begin
        for (shape, h) in (("linear", 0.3), ("exponential", 0.2))
            δ = 1e-6
            @test (radial_signal(R=R, D=D, times=[δ], layer=(shape, ρ, h))[1] - 1) / δ ≈ -2ρ / R rtol=1e-3
            @test radial_signal(R=R, D=D, times=[50.0], layer=(shape, ρ, h))[1] >= exp(-2ρ / R * 50) * (1 - 1e-9)
        end
    end
    @testset "grid convergence" begin
        for (shape, h) in (("linear", 0.3), ("exponential", 0.2))
            a = radial_signal(R=R, D=D, times=[50.0], layer=(shape, ρ, h), N=1000)[1]
            b = radial_signal(R=R, D=D, times=[50.0], layer=(shape, ρ, h), N=2000)[1]
            @test abs(a - b) < 2e-8
        end
    end
end
