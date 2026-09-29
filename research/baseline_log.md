# Baseline test log

## 2026-09-29

- Julia: 1.12.7
- Code commit tested: d9c5d0e
- Command: julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["evolve"])'
- Result: PASS — 113 of 113 tests passed in 16.9 seconds.
- Dependency issues: None observed during this test.

- Command: julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["radio_frequency"])'
- Result: PASS — 69 of 69 tests passed in 8.4 seconds.
- Dependency issues: None observed during this test.

