# P7 Steane Code — Sample Write-Up Outline (NOT A FINAL SUBMISSION)

Export your completed 2–4 page report as `WRITEUP.pdf` in the repository root.

## Title / Team
Steane [[7,1,3]] Quantum Error Correction — [TEAM NAME]

## 1. Problem
Single physical qubits can suffer X, Y and Z errors. The goal is to encode a logical qubit across seven physical qubits and detect/correct any one physical Pauli error.

## 2. Method
**Current sample:** classical parity-check syndrome algebra using a 3×7 Hamming matrix, with separate X and Z error-component syndromes. No circuit is built or executed.

**Replace with real implementation:** describe encoding circuit, six stabilizer generators, ancilla extraction, measurement/correction procedure, logical gate, simulator/noise parameters, number of shots and random seeds.

## 3. Results
**Illustrative starter check:** all 21 X/Y/Z single-position error *labels* map to themselves through algebraic syndrome lookup. This is not evidence of recovered logical qubit fidelity.

**TO COMPLETE:** include one measured/simulated figure, error bars if statistical, circuit-level verification against ideal reference, and a quantitative interpretation.

## 4. Honest Limits
Starter uses perfect parity checks and an assumed single error. It does not simulate quantum noise, noisy syndrome extraction, measurement faults or fault tolerance. Multi-error failure and thresholds remain untested.

## 5. AI use
See `AI_USE_NOTE.md` and replace placeholder details with verified, truthful disclosures.

## References
[Add exact papers, Qiskit docs and any borrowed code or algorithms.]
