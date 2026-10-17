"""Circuit-level experiments for the Steane code (run with Qiskit Aer)."""
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.quantum_info import Statevector, DensityMatrix, partial_trace, state_fidelity, Operator
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
import steane as S


# ----------------------------------------------------------------------------- Core task 3: 21 single-qubit errors
def correction_circuit(theta, phi, error=None):
    """Prepare cos(θ/2)|0>+e^{iφ}sin(θ/2)|1>, encode, inject `error`=(pauli, qubit), measure the six
    stabilisers with ancillas, and correct with mid-circuit feed-forward (if_test)."""
    data = QuantumRegister(7, "d"); anc = QuantumRegister(6, "a")
    sx = ClassicalRegister(3, "sx"); sz = ClassicalRegister(3, "sz")
    qc = QuantumCircuit(data, anc, sx, sz)
    qc.u(theta, phi, 0, data[S.ENCODER_INPUT])
    qc.compose(S.encoder_circuit(), data, inplace=True)
    qc.barrier()
    if error is not None:
        errs = [error] if isinstance(error, tuple) else error     # one error or a list of errors
        for pa, qq in errs:
            getattr(qc, pa.lower())(data[qq])
    qc.barrier()
    S.add_syndrome_extraction(qc, data, anc, sx, sz)
    S.add_feedforward_correction(qc, data, sx, sz)
    qc.save_statevector()
    return qc


def run_correction(theta, phi, error=None, seed=1):
    qc = correction_circuit(theta, phi, error)
    sim = AerSimulator(method="statevector", seed_simulator=seed)
    res = sim.run(transpile(qc, sim), shots=1).result()
    sv = Statevector(res.get_statevector())
    counts = res.get_counts()
    bits = list(counts.keys())[0].split()           # "sz sx" (last-added register first)
    sz_v, sx_v = int(bits[0], 2), int(bits[1], 2)
    rho_data = partial_trace(sv, list(range(7, 13)))  # trace out the ancillas
    ideal = Statevector(S.isometry() @ np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)]))
    fid = state_fidelity(rho_data, ideal)
    return fid, sx_v, sz_v


def single_error_table(test_states=None):
    """Rows for all 21 single-qubit errors (+ the no-error control)."""
    if test_states is None:
        rng = np.random.default_rng(7)
        test_states = [(0.0, 0.0), (np.pi, 0.0), (np.pi / 2, 0.0), (np.pi / 2, np.pi / 2),
                       (float(rng.uniform(0, np.pi)), float(rng.uniform(0, 2 * np.pi)))]
    rows = []
    for err in [None] + [(p, q) for q in range(7) for p in "XYZ"]:
        fids = []; syn = None
        for th, ph in test_states:
            f, sx_v, sz_v = run_correction(th, ph, err)
            fids.append(f); syn = (sx_v, sz_v)
        rows.append(dict(error="none" if err is None else f"{err[0]}{err[1]}",
                         sx=syn[0], sz=syn[1], min_fid=min(fids)))
    return rows


# ----------------------------------------------------------------------------- Core task 4: transversal gates
def logical_matrix_single(gate_name):
    """2x2 logical action of `gate^{⊗7}`: V† (U^{⊗7}) V, and the leakage out of the code space."""
    V = S.isometry()
    qc = QuantumCircuit(7)
    for q in range(7):
        getattr(qc, gate_name)(q)
    W = Operator(qc).data @ V
    M = V.conj().T @ W
    leak = np.linalg.norm(W - V @ M)
    return M, leak


def logical_cnot():
    """Logical 4x4 action of the transversal CNOT between two blocks (block A = control)."""
    V = S.isometry()
    qc = QuantumCircuit(14)
    for q in range(7):
        qc.cx(q, 7 + q)
    cols = []; leak = 0.0
    VV = lambda a, b: np.kron(V[:, b], V[:, a])      # qubits 0..6 = block A (low bits), 7..13 = block B
    for b in (0, 1):
        for a in (0, 1):
            out = Statevector(VV(a, b)).evolve(qc).data
            col = np.array([np.vdot(VV(a2, b2), out) for b2 in (0, 1) for a2 in (0, 1)])
            leak += np.linalg.norm(out - sum(col[i] * VV(i % 2, i // 2) for i in range(4)))
            cols.append(col)
    return np.array(cols).T, leak


def logical_bell_checks():
    """Encode |+>_L ⊗ |0>_L, apply the transversal CNOT, and check all 12 block stabilisers plus
    logical XX and ZZ: a logical Bell pair."""
    V = S.isometry()
    plus = V @ np.array([1, 1]) / np.sqrt(2); zero = V[:, 0]
    qc = QuantumCircuit(14)
    for q in range(7):
        qc.cx(q, 7 + q)
    sv = Statevector(np.kron(zero, plus)).evolve(qc)          # A = |+>_L (low bits), B = |0>_L
    out = {}
    for blk, off in (("A", 0), ("B", 7)):
        vals = []
        for name, P in S.stabilizers():
            lab = ["I"] * 14
            for q in range(7):
                ch = str(P)[::-1][q]
                lab[13 - (q + off)] = ch
            vals.append(sv.expectation_value(Pauli_from(lab)).real)
        out[blk] = vals
    xx = ["I"] * 14; zz = ["I"] * 14
    for q in range(14):
        xx[13 - q] = "X"; zz[13 - q] = "Z"
    out["XX_L"] = sv.expectation_value(Pauli_from(xx)).real
    out["ZZ_L"] = sv.expectation_value(Pauli_from(zz)).real
    return out


def Pauli_from(label_list):
    from qiskit.quantum_info import Pauli
    return Pauli("".join(label_list))


# ----------------------------------------------------------------------------- Stretch: noisy memory experiment
def make_noise_model(p, mode="circuit"):
    """mode='circuit': every gate and measurement is faulty (depolarising p on gates, readout flip p).
    mode='capacity': only the explicit `id` gate on the data qubits is faulty (matches the exact analysis)."""
    nm = NoiseModel()
    # NOTE: Aer's depolarizing_error(lam, n) is rho -> (1-lam) rho + lam I/2^n, so the probability of a
    # non-identity Pauli is lam*(1 - 4^-n), NOT lam. We rescale so that p is the total Pauli-error probability.
    lam1, lam2 = 4 * p / 3, 16 * p / 15
    if mode == "capacity":
        nm.add_all_qubit_quantum_error(depolarizing_error(lam1, 1), ["id"])
    else:
        nm.add_all_qubit_quantum_error(depolarizing_error(lam1, 1), ["h", "id"])
        nm.add_all_qubit_quantum_error(depolarizing_error(lam2, 2), ["cx"])
        nm.add_all_qubit_readout_error(ReadoutError([[1 - p, p], [p, 1 - p]]))
    return nm


def memory_circuit(basis, rounds, mode="circuit"):
    """Logical memory: encode |0>_L (basis 'Z') or |+>_L (basis 'X'), `rounds` rounds of syndrome
    extraction with fresh ancillas, then transversal readout of all 7 data qubits.
    Classical bits: 0..6 data readout; 7+6r.. : round r (X-type rows 0..2, Z-type rows 3..5)."""
    nq = 7 + 6 * rounds
    qc = QuantumCircuit(nq, 7 + 6 * rounds)
    if basis == "X":
        qc.h(S.ENCODER_INPUT)
    qc.compose(S.encoder_circuit(), range(7), inplace=True)
    if mode == "capacity":
        for q in range(7):
            qc.id(q)                                  # the only noisy location
    for r in range(rounds):
        anc = [7 + 6 * r + i for i in range(6)]
        base = 7 + 6 * r
        for i in range(3):
            qc.h(anc[i])
            for q in S.SUPPORTS[i]:
                qc.cx(anc[i], q)
            qc.h(anc[i]); qc.measure(anc[i], base + i)
        for i in range(3):
            for q in S.SUPPORTS[i]:
                qc.cx(q, anc[3 + i])
            qc.measure(anc[3 + i], base + 3 + i)
    if basis == "X":
        for q in range(7):
            qc.h(q)
    for q in range(7):
        qc.measure(q, q)
    return qc


def _val(bits3):
    return int(bits3[0] * 4 + bits3[1] * 2 + bits3[2])


def decode_shots(mem, basis, rounds, decoder="frame"):
    """Classical decoding of the shots. decoder='frame': apply each round's correction to a Pauli frame
    (identical in effect to in-circuit feed-forward for Pauli corrections) then Hamming-decode the final
    readout. decoder='agree': only trust a syndrome that repeats in two consecutive rounds."""
    fails = 0
    for s in mem:
        b = np.array([int(c) for c in s[::-1]])           # b[i] = classical bit i
        m = b[:7].copy()
        off = 3 if basis == "Z" else 0                     # Z-basis flips come from Z-type stabilisers (sz)
        synd = [_val(b[7 + 6 * r + off: 7 + 6 * r + off + 3]) for r in range(rounds)]
        if decoder == "agree":
            for r in range(1, rounds):
                if synd[r] and synd[r] == synd[r - 1]:
                    m[synd[r] - 1] ^= 1
                    synd = synd[:r] + [0] * (rounds - r)   # consumed
        else:
            for v in synd:
                if v:
                    m[v - 1] ^= 1
        fv = _val(S.H_MAT.dot(m) % 2)                      # final classical Hamming decode of the readout
        if fv:
            m[fv - 1] ^= 1
        fails += int(m.sum() % 2)                          # logical value (0 expected in both bases)
    return fails


def memory_failure_rate(p, rounds, mode="circuit", shots=20000, seed=11, decoder="frame"):
    """Average of the Z-basis and X-basis logical failure probability, plus its binomial standard error."""
    sim = AerSimulator(method="stabilizer", noise_model=make_noise_model(p, mode), seed_simulator=seed)
    rates = []
    for basis in ("Z", "X"):
        qc = transpile(memory_circuit(basis, rounds, mode), sim, optimization_level=0)
        res = sim.run(qc, shots=shots, memory=True).result()
        f = decode_shots(res.get_memory(), basis, rounds, decoder)
        rates.append(f / shots)
    mean = float(np.mean(rates)); err = float(np.sqrt(sum(r * (1 - r) / shots for r in rates)) / 2)
    return mean, err, rates


def bare_failure_rate(p, shots=20000, seed=5):
    """Baseline: one physical qubit, one noisy idle gate, noisy readout (same noise model, Z and X basis)."""
    sim = AerSimulator(method="stabilizer", noise_model=make_noise_model(p, "circuit"), seed_simulator=seed)
    rates = []
    for basis in ("Z", "X"):
        qc = QuantumCircuit(1, 1)
        if basis == "X": qc.h(0)
        qc.id(0)
        if basis == "X": qc.h(0)
        qc.measure(0, 0)
        res = sim.run(transpile(qc, sim, optimization_level=0), shots=shots).result().get_counts()
        rates.append(res.get("1", 0) / shots)
    mean = float(np.mean(rates)); err = float(np.sqrt(sum(r * (1 - r) / shots for r in rates)) / 2)
    return mean, err
