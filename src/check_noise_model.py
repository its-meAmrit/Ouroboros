"""Independent check that Aer's noise sampling reproduces the exact analysis (no code circuit involved).
Applies noisy identity gates to 7 qubits, decodes the resulting bit flips with the Hamming lookup, and compares the
failure rate with the exact enumeration. Usage: python src/check_noise_model.py   (about 1 minute)"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import experiments as E, steane as S
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

p, shots = 0.1, 1_000_000
sim = AerSimulator(method="stabilizer", noise_model=E.make_noise_model(p, "capacity"), seed_simulator=123)
qc = QuantumCircuit(7, 7)
for q in range(7): qc.id(q)
qc.measure(range(7), range(7))
mem = sim.run(transpile(qc, sim, optimization_level=0), shots=shots, memory=True).result().get_memory()
B = np.array([[int(c) for c in s[::-1]] for s in mem], dtype=np.int8)
s = (B @ S.H_MAT.T) % 2; v = s[:, 0] * 4 + s[:, 1] * 2 + s[:, 2]
idx = np.flatnonzero(v); B[idx, v[idx] - 1] ^= 1
f = (B.sum(1) % 2).mean()
print(f"Aer: {f:.5f} +- {(f * (1 - f) / shots) ** .5:.5f}    exact: {S.exact_logical_error(p, 'x'):.5f}")
