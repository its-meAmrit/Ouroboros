"""Steane [[7,1,3]] code in Qiskit: stabilisers, encoder, syndrome extraction,
lookup-table decoder, transversal gates, and exact / circuit-level error analysis.

Qubit convention: data qubits are numbered 0..6. Position j (0-indexed) corresponds
to column j+1 of the Hamming parity-check matrix, i.e. to the 3-bit binary number j+1.
That choice makes the syndrome, read as a binary number, equal to (error position + 1).
"""
import itertools
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.quantum_info import Pauli, Operator, Statevector

N = 7
# Hamming parity-check matrix: row r, column j is bit (2-r) of the number (j+1).
H_MAT = np.array([[((j + 1) >> (2 - r)) & 1 for j in range(N)] for r in range(3)], dtype=int)
# -> [[0 0 0 1 1 1 1],
#     [0 1 1 0 0 1 1],
#     [1 0 1 0 1 0 1]]
SUPPORTS = [tuple(int(j) for j in np.flatnonzero(H_MAT[r])) for r in range(3)]

ENCODER_INPUT = 2                         # qubit that carries the unencoded state
PIVOTS = {3: (4, 5, 6), 1: (2, 5, 6), 0: (2, 4, 6)}   # pivot -> targets (one per X stabiliser)


def pauli_on(kind, support):
    """Pauli string with `kind` ('X' or 'Z') on the given qubits (Qiskit label is little-endian)."""
    label = ["I"] * N
    for q in support:
        label[N - 1 - q] = kind
    return Pauli("".join(label))


def stabilizers():
    """The six stabiliser generators as (name, Pauli). X-type detect Z errors, Z-type detect X errors."""
    out = []
    for r in range(3):
        out.append((f"gX{r+1}", pauli_on("X", SUPPORTS[r])))
    for r in range(3):
        out.append((f"gZ{r+1}", pauli_on("Z", SUPPORTS[r])))
    return out


LOGICAL_X = pauli_on("X", range(N))      # X^{⊗7}
LOGICAL_Z = pauli_on("Z", range(N))      # Z^{⊗7}


# --------------------------------------------------------------------------- encoder
def encoder_circuit():
    """Encode the state on qubit 2 into the 7-qubit Steane block (other qubits start in |0>)."""
    qc = QuantumCircuit(N, name="Steane encoder")
    qc.cx(2, 4)                           # a|0>+b|1> -> a|000>+b|111> on {2,4,5}; {2,4,5} is a weight-3 codeword,
    qc.cx(2, 5)                           # so X_{245} is a logical X (equivalent to X^{⊗7} on the code space)
    for p in PIVOTS:
        qc.h(p)                           # pivots into |+>: each will control one X-stabiliser
    for p, targets in PIVOTS.items():
        for t in targets:
            qc.cx(p, t)                   # controlled application of (1 + g_X) -> projects onto X-stabiliser +1
    return qc


def isometry():
    """128x2 matrix V with V|0> = |0_L>, V|1> = |1_L> (built from the encoder circuit)."""
    U = Operator(encoder_circuit()).data
    return U[:, [0, 1 << ENCODER_INPUT]]


# --------------------------------------------------------------------------- syndrome + correction
def add_syndrome_extraction(qc, data, anc, sx, sz):
    """Measure all six stabilisers with one ancilla each.
    anc[0..2]: X-type (detect Z errors) -> sx; anc[3..5]: Z-type (detect X errors) -> sz.
    Row r is stored in classical bit (2-r) so the register's integer value = error position + 1."""
    for r in range(3):
        a = anc[r]
        qc.h(a)
        for q in SUPPORTS[r]:
            qc.cx(a, data[q])
        qc.h(a)
        qc.measure(a, sx[2 - r])
    for r in range(3):
        a = anc[3 + r]
        for q in SUPPORTS[r]:
            qc.cx(data[q], a)
        qc.measure(a, sz[2 - r])


def add_feedforward_correction(qc, data, sx, sz):
    """Lookup table, executed inside the circuit with if_test (mid-circuit feed-forward)."""
    for v in range(1, 8):
        with qc.if_test((sz, v)):
            qc.x(data[v - 1])             # Z-stabiliser syndrome v  -> X error on qubit v-1
        with qc.if_test((sx, v)):
            qc.z(data[v - 1])             # X-stabiliser syndrome v  -> Z error on qubit v-1


def correction_table():
    """Human-readable lookup table: syndrome value -> position."""
    return {v: v - 1 for v in range(1, 8)}


# --------------------------------------------------------------------------- exact (code-capacity) analysis
def _syndrome_val(vec):
    s = H_MAT.dot(vec) % 2
    return int(s[0] * 4 + s[1] * 2 + s[2])


def residual(ex, ez):
    """Apply the lookup-table decoder to the error (ex, ez) and return the residual (rx, rz)."""
    rx, rz = np.array(ex), np.array(ez)
    sz, sx = _syndrome_val(ex), _syndrome_val(ez)
    if sz:
        rx = rx.copy(); rx[sz - 1] ^= 1
    if sx:
        rz = rz.copy(); rz[sx - 1] ^= 1
    return rx, rz


def is_logical_failure(ex, ez):
    """Residual lies in the normaliser by construction, so it is a stabiliser (harmless) iff both
    parts have even weight; odd weight means a logical X and/or Z has been applied."""
    rx, rz = residual(ex, ez)
    return bool(rx.sum() % 2), bool(rz.sum() % 2)       # (logical-X flip, logical-Z flip)


def all_paulis():
    """Yield every one of the 4^7 Pauli errors as (ex, ez, weight)."""
    for bits in itertools.product(range(4), repeat=N):   # 0=I 1=X 2=Z 3=Y
        ex = np.array([b & 1 for b in bits]); ez = np.array([(b >> 1) & 1 for b in bits])
        yield ex, ez, int(np.count_nonzero(bits))


_ENUM = None
def _enumeration():
    global _ENUM
    if _ENUM is None:
        rows = []
        for ex, ez, w in all_paulis():
            fx, fz = is_logical_failure(ex, ez)
            rows.append((w, fx, fz))
        _ENUM = np.array(rows, dtype=int)
    return _ENUM


def failing_counts():
    """number of Pauli errors of each weight w that cause a logical error: dict w -> (failing, total)."""
    E = _enumeration(); out = {}
    for w in range(N + 1):
        sel = E[E[:, 0] == w]
        out[w] = (int(np.count_nonzero(sel[:, 1] | sel[:, 2])), len(sel))
    return out


def exact_logical_error(p, which="any"):
    """Exact logical error probability under iid depolarising noise (X,Y,Z each p/3 per qubit),
    perfect syndrome extraction and lookup-table decoding. which: 'any' | 'x' | 'z'."""
    E = _enumeration()
    w = E[:, 0]
    prob = (p / 3) ** w * (1 - p) ** (N - w)
    if which == "x":
        mask = E[:, 1] == 1
    elif which == "z":
        mask = E[:, 2] == 1
    else:
        mask = (E[:, 1] | E[:, 2]) == 1
    return float(prob[mask].sum())
