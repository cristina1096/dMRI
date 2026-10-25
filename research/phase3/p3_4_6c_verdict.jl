# P3.4.6 verdict over both Monte Carlo passes (comparison.csv: τ_a, τ_b; comparison_tau16.csv: τ_c for conventional,
# h = 0.05, 0.1). Per model the criterion (max over readouts |z| ≤ 3 for S and for the paired ΔS vs the continuum) is
# applied at the smallest τ run for it. h = 1.0 is deterministic in MC (uniform rate, SEM ≈ 1e-17), so its S is compared
# exactly, |S_MC/S_cont − 1| ≤ 2e-6 (continuum roundoff bound), instead of a z. For the paired ΔS the timestep trend is
# reported with a Richardson estimate ΔS(τ → 0) ≈ ΔS_c − (ΔS_b − ΔS_c)/3 (linear in τ; τ_b/τ_c = 4).
# Run: julia --project=research/baseline research/phase3/p3_4_6c_verdict.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_6_continuum_check")
readrows(f) = (ls = readlines(joinpath(OUT, f)); h = split(ls[1], ','); [Dict(zip(h, split(l, ','))) for l in ls[2:end]])
num(r, k) = parse(Float64, r[k])
first, third = readrows("comparison.csv"), readrows("comparison_tau16.csv")
τb = minimum(num(r, "tau") for r in first)
verdict = Dict{String, Any}()
for r in filter(r -> num(r, "tau") == τb, first)
    name = r["model"]
    c = filter(x -> x["model"] == name, third)
    if !isempty(c)
        q = only(c)
        rich = num(q, "dS50_mc") - (num(r, "dS50_mc") - num(q, "dS50_mc")) / 3
        verdict[name] = Dict("tau" => num(q, "tau"), "maxabs_zS" => num(q, "maxabs_zS"), "maxabs_zdS" => num(q, "maxabs_zdS"),
            "dS50_tau_b" => num(r, "dS50_mc"), "dS50_tau_c" => num(q, "dS50_mc"), "dS50_richardson" => rich,
            "dS50_cont" => num(q, "dS50_cont"), "pass" => num(q, "maxabs_zS") <= 3 && num(q, "maxabs_zdS") <= 3)
    elseif name == "step_h1.0"
        rel = abs(num(r, "S50_mc") / num(r, "S50_cont") - 1)
        verdict[name] = Dict("tau" => τb, "S50_rel_dev" => rel, "maxabs_zdS" => num(r, "maxabs_zdS"),
            "pass" => rel <= 2e-6 && num(r, "maxabs_zdS") <= 3)
    else
        verdict[name] = Dict("tau" => τb, "maxabs_zS" => num(r, "maxabs_zS"), "maxabs_zdS" => num(r, "maxabs_zdS"),
            "pass" => num(r, "maxabs_zS") <= 3 && num(r, "maxabs_zdS") <= 3)
    end
end
for k in sort(collect(keys(verdict)))
    println(rpad(k, 14), verdict[k])
end
allpass = all(v["pass"] for v in values(verdict))
println(allpass ? "PASS" : "FAIL (" * join([k for k in keys(verdict) if !verdict[k]["pass"]], ", ") * ")")
write_json(joinpath(OUT, "summary_final.json"), Dict("provenance" => provenance(), "pass" => allpass, "verdict" => verdict))
