"""Checks against exact references. Run:  pytest -q"""
import os, sys
import numpy as np
from qiskit.quantum_info import Statevector
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import steane as S, experiments as E


def test_stabilisers_commute_and_logicals_anticommute():
    gens = [P for _, P in S.stabilizers()]
    assert all(a.commutes(b) for a in gens for b in gens)
    assert all(g.commutes(S.LOGICAL_X) and g.commutes(S.LOGICAL_Z) for g in gens)
    assert not S.LOGICAL_X.commutes(S.LOGICAL_Z)


def test_encoder_gives_codewords():
    V = S.isometry()
    for b in (0, 1):
        for _, P in S.stabilizers():
            assert abs(Statevector(V[:, b]).expectation_value(P) - 1) < 1e-9
    # |0_L> is an equal superposition of the 8 even-weight Hamming codewords
    amp = np.abs(V[:, 0]); nz = np.flatnonzero(amp > 1e-9)
    assert len(nz) == 8 and np.allclose(amp[nz], 1 / np.sqrt(8))
    assert all(bin(i).count("1") % 2 == 0 and not (S.H_MAT @ [(i >> q) & 1 for q in range(7)] % 2).any() for i in nz)


def test_all_21_single_errors_corrected_by_circuit():
    rows = E.single_error_table(test_states=[(0.3, 1.0), (np.pi / 2, 0.0)])
    assert len(rows) == 22 and all(r["min_fid"] > 1 - 1e-9 for r in rows)
    syn = {(r["sx"], r["sz"]) for r in rows[1:]}
    assert len(syn) == 21                         # every single-qubit error has its own syndrome pair


def test_two_errors_fail():
    assert E.run_correction(1.1, 0.4, [("X", 0), ("X", 1)])[0] < 0.99
    assert E.run_correction(1.1, 0.4, [("X", 0), ("Z", 1)])[0] > 1 - 1e-9     # different types are fine


def test_failing_weight_counts():
    fc = S.failing_counts()
    assert fc[1][0] == 0 and fc[2] == (147, 189)


def test_transversal_gates():
    s2 = 1 / np.sqrt(2)
    for g, M in {"x": [[0, 1], [1, 0]], "z": [[1, 0], [0, -1]], "h": [[s2, s2], [s2, -s2]],
                 "sdg": [[1, 0], [0, 1j]]}.items():
        got, leak = E.logical_matrix_single(g)
        assert np.allclose(got, M) and leak < 1e-9, g
    got, leak = E.logical_matrix_single("s")            # S^{⊗7} is logical S^dagger
    assert np.allclose(got, [[1, 0], [0, -1j]]) and leak < 1e-9
    M, leak = E.logical_cnot()
    assert np.allclose(M, [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]) and leak < 1e-9


def test_exact_matches_circuit_sim_and_limits():
    assert S.exact_logical_error(0.0) == 0
    p = 0.1
    m, e, _ = E.memory_failure_rate(p, 1, "capacity", shots=30000, seed=3)
    assert abs(m - S.exact_logical_error(p, "x")) < 4 * e + 0.003          # agree within statistics
    assert abs(S.exact_logical_error(0.01) - 147 / 9 * 0.01 ** 2) < 3e-4   # leading-order p^2 law
