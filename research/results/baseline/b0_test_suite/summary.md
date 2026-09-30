# B0 — existing test suite on unmodified code

- Command: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["collisions","evolve","permeability","transfer","known_sequences"])'`
- Date: 2026-09-30; branch `cc/near-surface-r2`, HEAD `677a708`, `src/` unmodified
- Julia 1.12.7, `julia` default threads, Apple M1 Pro
- Result: **442 / 442 pass**, 0 fail, 0 error, 0 broken (test time 4 min 36 s; wall clock 9 min 23 s incl. precompilation)

| Suite | Time (s) |
|---|---|
| collisions | 55.3 |
| evolve | 5.7 |
| permeability | 24.8 |
| transfer | 54.9 |
| known_sequences | 135.6 |

Full output: `run.log`.
