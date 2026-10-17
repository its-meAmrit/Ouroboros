# AI-use note

## Tools used

Claude (Anthropic, claude.ai chat with a code sandbox) wrote the first version of all code in `src/`, the tests, and the write-up, and ran the simulations. Qiskit/Aer documentation was not consulted directly; behaviour was checked by experiment instead.

## One thing the AI got wrong

The first noise model passed `p` straight into Aer's `depolarizing_error(p, n)`. That parameter is the mixing weight λ in ρ → (1−λ)ρ + λ·I/2ⁿ, so the probability of an actual non-identity Pauli is λ(1 − 4⁻ⁿ), i.e. only 3p/4 for one qubit (15p/16 for two), not p. The code ran without any error and produced smooth, plausible curves, but they were silently shifted to lower noise.

## How we caught it

A cross-check of the Aer circuit simulation against our exact enumeration of all 4⁷ Pauli errors. At p = 0.1 the logical failure probability from Aer was 0.0406 per basis, while the exact value was 0.0681 (a mismatch of more than 25 standard errors). Plugging an effective p = 0.075 into the exact formula gives 0.0415, which pinned the discrepancy on the 3/4 factor. After rescaling (λ = 4p/3 and 16p/15) the two agree. A separate 1,000,000-shot run with no code circuit at all (`src/check_noise_model.py`) gives 0.06824 ± 0.00025 against an exact 0.06813. A second, smaller slip: the first version of Fig. (a) compared one-basis Aer points with an any-logical-error curve, so the points looked too low; this was spotted by looking at the plot and fixed by adding the matching one-basis exact curve.

## What we verified ourselves

The `tests/` suite checks every circuit against an exact reference: stabilisers commute, the encoder output is the uniform superposition of the 8 even-weight Hamming codewords with all six stabilisers at +1, all 21 single errors are corrected with fidelity 1 on five states, 147/189 weight-2 errors fail (counted exhaustively), transversal gates reproduce the expected 2×2 and 4×4 logical matrices with zero leakage out of the code space, and the circuit simulation matches the exact curve within statistics. The write-up's Method section explains each circuit, formula and plot.
