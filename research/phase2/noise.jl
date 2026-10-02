"""
    required_spins(nspins, separation_in_sem; target=10)

Spins per seed needed so that a separation currently equal to `separation_in_sem` SEM (measured with
`nspins` spins per seed) becomes at least `target` SEM. SEM ∝ 1/√N, so N_new = N·(target/separation)².
Never recommends fewer spins than `nspins`.
"""
function required_spins(nspins::Integer, separation_in_sem::Real; target=10)
    separation_in_sem >= target && return Int(nspins)
    return ceil(Int, nspins * (target / separation_in_sem)^2)
end
