# P2.5.9 — overlapping layers (linear, exponential) vs the conventional model at the same ρ, exact continuum benchmark.
# (1) Overlap grid × ρw/D: max |ΔS| over the window (0 until the conventional signal ≈ 1 %, 41 readouts, R2_bulk = 0),
#     ΔT2/T2 at TE* = the conventional model's own T2 (S_conv(TE*) = e^−1; R2_bulk = 1/2000; ruling: a fixed TE = 50 is
#     unmeasurable at strong ρ) with the two-point noise at SNR 100/200, best-fit conventional residual
#     (ρ_c free, scan e^−4…e^8 × ρ + golden; single condition), and the layer's mono-exponential residual.
# (2) D dependence: apparent T2 (TE = TE* of the conventional model at D = 3, one per ρ) vs D ∈ {0.5, 1, 2, 3} at ρ = 0.1 and 1 for the layer at each overlap and for
#     the conventional model; relative change T2(D = 0.5)/T2(D = 3) − 1 per model.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_overlap.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))
include(joinpath(@__DIR__, "continuum.jl"))
include(joinpath(@__DIR__, "overlap.jl"))

const OUT = phase2_outdir("p2_5_overlap")
const W, D0 = 2.0, 3.0
const HW = Dict("linear" => [0.25, 0.5, 0.6, 0.75, 0.9, 1.0], "exponential" => [0.05, 0.1, 0.25, 0.5, 0.75, 1.0])
const RWD = [0.0067, 0.1, 0.3, 1.0, 3.0]

function best_fit(Sl, w, D, ρ0, times)
    f(lr) = maximum(abs.(Sl .- continuum_signal(w=w, D=D, times=times, surface=exp(lr), N=1500)))
    grid = log(ρ0) .+ (-4.0:0.25:8.0)
    v = f.(grid); k = argmin(v)
    a, b = grid[max(k - 1, 1)], grid[min(k + 1, length(grid))]
    gr = (sqrt(5) - 1) / 2
    x1, x2 = b - gr * (b - a), a + gr * (b - a)
    f1, f2 = f(x1), f(x2)
    for _ in 1:30
        if f1 < f2
            b, x2, f2 = x2, x1, f1
            x1 = b - gr * (b - a); f1 = f(x1)
        else
            a, x1, f1 = x1, x2, f2
            x2 = a + gr * (b - a); f2 = f(x2)
        end
    end
    return (resid=min(f1, f2, v[k]), at_edge=k == 1 || k == length(grid))
end
"TE* where the conventional signal (bulk 1/2000) equals e^(−1): its own apparent T2. Used as the echo time so that S(TE) stays measurable at every ρ."
function te_star(w, D, ρ)
    f(t) = continuum_signal(w=w, D=D, times=[t], surface=ρ, R2_bulk=1 / 2000)[1] - exp(-1)
    lo, hi = 1e-3, 1.0
    while f(hi) > 0
        hi *= 2
    end
    for _ in 1:50
        mid = (lo + hi) / 2
        f(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
const TESTAR = Dict{Tuple{Float64, Float64, Float64}, Float64}()
testar(w, D, ρ) = get!(() -> te_star(w, D, ρ), TESTAR, (w, D, ρ))

monoexp(t, S) = (k = findall(>(0), t); R = -sum(t[k] .* log.(S[k])) / sum(t[k] .^ 2); maximum(abs.(S .- exp.(-R .* t))))

rows = NamedTuple[]
for shape in ("linear", "exponential"), hw in HW[shape], rwd in RWD
    ρ, h = rwd * D0 / W, hw * W
    tmax = log(100) / continuum_rates(w=W, D=D0, surface=ρ, N=1500)[1]
    t = collect(range(0.0, tmax, length=41))
    Sl = continuum_signal(w=W, D=D0, times=t, layer=(shape, ρ, h))
    Sc = continuum_signal(w=W, D=D0, times=t, surface=ρ)
    bf = best_fit(continuum_signal(w=W, D=D0, times=t, layer=(shape, ρ, h), N=1500), W, D0, ρ, t)
    TE = testar(W, D0, ρ)
    St = [continuum_signal(w=W, D=D0, times=[0.0, TE], layer=(shape, ρ, h), R2_bulk=1 / 2000), continuum_signal(w=W, D=D0, times=[0.0, TE], surface=ρ, R2_bulk=1 / 2000)]
    T2l, T2c = apparent_t2(St[1][1], St[1][2], TE), apparent_t2(St[2][1], St[2][2], TE)
    r = St[2][2]
    m = overlap_metrics(shape, h, W)
    push!(rows, (shape=shape, h_over_w=hw, rho_w_over_D=rwd, rho=ρ, overlap_width_um=m.overlap_width, mid_over_wall=m.mid_over_wall,
        centroid_um=centroid(shape, h; cutoff=W), max_abs_dS=maximum(abs.(Sl .- Sc)), dS_sign=sign(Sl[argmax(abs.(Sl .- Sc))] - Sc[argmax(abs.(Sl .- Sc))]),
        TE_ms=TE, rel_T2=(T2l - T2c) / T2c, two_sigma_T2_snr200=2t2_noise_rel(r, 200), two_sigma_T2_snr100=2t2_noise_rel(r, 100),
        best_fit_residual=bf.resid, best_fit_at_edge=bf.at_edge, layer_monoexp_residual=monoexp(t, Sl), conv_monoexp_residual=monoexp(t, Sc)))
    q = rows[end]
    @printf("%-11s h/w=%.2f ρw/D=%-6g: mid/wall %.2e, max|ΔS| %.2e, ΔT2/T2 %+.2e at TE* = %.3g ms (2σ200 %.3f), best-fit %.2e%s, mono-exp layer %.1e conv %.1e\n",
        shape, hw, rwd, q.mid_over_wall, q.max_abs_dS, q.rel_T2, q.TE_ms, q.two_sigma_T2_snr200, q.best_fit_residual, q.best_fit_at_edge ? " (edge)" : "", q.layer_monoexp_residual, q.conv_monoexp_residual)
end
write_csv(joinpath(OUT, "grid.csv"), rows)

ddep = NamedTuple[]
for ρ in (0.1, 1.0)
    TE = testar(W, D0, ρ)                                       # one echo time per ρ (conventional T2 at D = 3), same for every D
    T2(model, D) = (S = model == "conventional" ? continuum_signal(w=W, D=D, times=[0.0, TE], surface=ρ, R2_bulk=1 / 2000) :
                        continuum_signal(w=W, D=D, times=[0.0, TE], layer=(model[1], ρ, model[2] * W), R2_bulk=1 / 2000); apparent_t2(S[1], S[2], TE))
    models = vcat(["conventional"], [(s, hw) for s in ("linear", "exponential") for hw in HW[s]])
    for mdl in models
        Ts = [T2(mdl, D) for D in (0.5, 1.0, 2.0, 3.0)]
        push!(ddep, (rho=ρ, TE_ms=TE, model=mdl == "conventional" ? "conventional" : mdl[1], h_over_w=mdl == "conventional" ? 0.0 : mdl[2],
            T2_D05=Ts[1], T2_D1=Ts[2], T2_D2=Ts[3], T2_D3=Ts[4], rel_change_D05_vs_D3=Ts[1] / Ts[4] - 1))
        @printf("D-dependence ρ=%.1f %-12s h/w=%-5s: T2(D=0.5..3) = %s ms; T2(0.5)/T2(3) − 1 = %+.2e\n", ρ, ddep[end].model,
            mdl == "conventional" ? "-" : string(mdl[2]), join(round.(Ts, digits=3), ", "), ddep[end].rel_change_D05_vs_D3)
    end
end
write_csv(joinpath(OUT, "ddep.csv"), ddep)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "w" => W, "D" => D0, "rows" => rows, "ddep" => ddep))
println("saved to $OUT")
