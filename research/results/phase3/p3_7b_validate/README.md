# p3_7b_validate (raw Monte Carlo output)

The `S_sup` / `dS_sup` columns of `curves.csv` hold the annulus continuum **as it was at run time: cut at R_out and NOT
renormalised**. After the validation the annulus layer is renormalised by 1/kept by default (user decision 2026-11-09,
`continuum_annulus.jl`); `p3_7b_validate_analysis.jl` recomputes the renormalised annulus and keeps the run's values as
`*_cut` for transparency. Later tasks (inside, density, map, ddep) write the renormalised annulus in `S_sup`.
