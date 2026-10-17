"""Reproduces every number and figure of the Steane-code project (Problem 7, Track 4).
Usage:  python src/main.py [--quick]      (--quick: fewer shots, ~1 min; default ~4 min)"""
import argparse, json, os, sys, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(__file__))
import steane as S
import experiments as E

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)
SEED = 42


def fig_encoder():
    qc = S.encoder_circuit()
    try:
        qc.draw("mpl", filename=os.path.join(FIG, "encoder_circuit.png"), fold=-1)
    except Exception as e:                                   # fall back to text rendering
        txt = str(qc.draw("text", fold=-1))
        fig = plt.figure(figsize=(9, 3.2)); fig.text(0.01, 0.5, txt, family="monospace", fontsize=6, va="center")
        fig.savefig(os.path.join(FIG, "encoder_circuit.png"), dpi=200); plt.close(fig)


def fig_table(rows):
    fig, ax = plt.subplots(figsize=(9, 5.8)); ax.axis("off")
    cell = [[r["error"], f"{r['sx']} ({r['sx']:03b})", f"{r['sz']} ({r['sz']:03b})",
             ("-" if r["error"] == "none" else
              "  ".join(filter(None, [f"X{r['sz']-1}" if r['sz'] else "", f"Z{r['sx']-1}" if r['sx'] else ""]))),
             f"{r['min_fid']:.12f}", "PASS" if r["min_fid"] > 1 - 1e-9 else "FAIL"] for r in rows]
    tb = ax.table(cellText=cell, colLabels=["injected", "X-stab. syndrome", "Z-stab. syndrome", "correction applied",
                                            "min fidelity (5 states)", ""], loc="center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 1.15)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "single_error_table.png"), dpi=200); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); a = ap.parse_args()
    shots = 4000 if a.quick else 10000
    t0 = time.time(); out = {"seed": SEED, "shots": shots}

    # ---- Core 1: stabilisers / encoder
    V = S.isometry()
    out["stabilisers"] = {n: str(P) for n, P in S.stabilizers()}
    from qiskit.quantum_info import Statevector
    out["stabiliser_expectations_on_0L_1L"] = {n: [float(Statevector(V[:, b]).expectation_value(P).real) for b in (0, 1)]
                                              for n, P in S.stabilizers()}
    print("encoder:", S.encoder_circuit().count_ops())
    fig_encoder()

    # ---- Core 2+3: syndrome extraction, lookup, all 21 single errors (with mid-circuit feed-forward)
    rows = E.single_error_table()
    out["single_error_table"] = rows
    out["all_21_corrected"] = all(r["min_fid"] > 1 - 1e-9 for r in rows)
    print("21 single errors all corrected:", out["all_21_corrected"])
    fig_table(rows)
    with open(os.path.join(ROOT, "figures", "single_error_table.csv"), "w") as f:
        f.write("error,x_syndrome,z_syndrome,min_fidelity\n")
        for r in rows: f.write(f"{r['error']},{r['sx']},{r['sz']},{r['min_fid']:.12f}\n")

    # ---- two errors fail
    two = {}
    for name, errs in {"X0 X1 (two bit flips)": [("X", 0), ("X", 1)], "Z2 Z5 (two phase flips)": [("Z", 2), ("Z", 5)],
                       "X0 Z1 (one of each, different qubits)": [("X", 0), ("Z", 1)],
                       "Y0 Y1": [("Y", 0), ("Y", 1)]}.items():
        two[name] = E.run_correction(1.1, 0.4, errs)[0]
    out["two_error_fidelity"] = two; print("two-error fidelities:", two)
    out["failing_pauli_counts_by_weight"] = {str(k): v for k, v in S.failing_counts().items()}

    # ---- Core 4: transversal gates
    tg = {}
    for g in ["x", "z", "h", "s", "sdg"]:
        M, leak = E.logical_matrix_single(g); tg[g] = {"logical_matrix": np.round(M, 6).tolist() and
                  [[complex(np.round(c, 6)).__repr__() for c in row] for row in M], "leakage": float(leak)}
    Mc, leakc = E.logical_cnot(); tg["cnot"] = {"logical_matrix": np.round(Mc.real, 6).tolist(), "leakage": float(leakc)}
    out["transversal"] = tg; out["logical_bell"] = E.logical_bell_checks()
    print("transversal CNOT logical matrix:\n", np.round(Mc.real, 3))

    # ---- Stretch: logical vs physical error rate
    ps = np.logspace(-3, np.log10(0.3), 60)
    exact = np.array([S.exact_logical_error(p) for p in ps])
    p_be = brentq(lambda x: S.exact_logical_error(x) - x, 0.02, 0.4)
    A = S.failing_counts()[2][0] / 9
    out["break_even_code_capacity"] = p_be; out["leading_coefficient_A"] = A
    print(f"code-capacity break-even p* = {p_be:.4f}; leading term p_L ~ {A:.2f} p^2")
    cap_p = [0.02, 0.05, 0.1, 0.15]
    cap = [E.memory_failure_rate(p, 1, "capacity", shots=shots * 2, seed=SEED) for p in cap_p]
    out["capacity_aer"] = [{"p": p, "mean": c[0], "stderr": c[1], "exact_xz_each": S.exact_logical_error(p, "x")}
                           for p, c in zip(cap_p, cap)]

    cps = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.04]
    circ = {}
    for label, R, dec in (("R=0 (encode + readout only)", 0, "frame"), ("R=1 round", 1, "frame"),
                          ("R=2 rounds, act only if syndromes agree", 2, "agree")):
        circ[label] = [E.memory_failure_rate(p, R, "circuit", shots=shots, seed=SEED, decoder=dec)[:2] for p in cps]
        print("done", label, round(time.time() - t0), "s")
    bare = [E.bare_failure_rate(p, shots=shots * 2, seed=SEED) for p in cps]
    out["circuit_level"] = {"p": cps, "bare": bare, **{k: v for k, v in circ.items()}}

    # ---- main figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    a0 = ax[0]
    a0.loglog(ps, exact, "C0-", lw=2, label="Steane, exact (all $4^7$ Paulis)")
    a0.loglog(ps, ps, "k--", label="unencoded qubit ($p_L=p$)")
    a0.loglog(ps[ps < 0.1], A * ps[ps < 0.1] ** 2, "C0:", label=f"leading term ${A:.1f}\\,p^2$")
    exact_x = np.array([S.exact_logical_error(p, "x") for p in ps])
    a0.loglog(ps, exact_x, "C3--", lw=1, label="Steane, exact, one basis only (what Aer measures)")
    a0.errorbar(cap_p, [c[0] for c in cap], yerr=[c[1] for c in cap], fmt="o", color="C3", capsize=3,
                label="Aer circuit simulation (Z/X-basis mean)")
    a0.plot([p_be], [p_be], "k*", ms=13); a0.annotate(f"break-even\n$p^*\\approx{p_be:.3f}$", (p_be, p_be), (p_be * 0.15, p_be * 3),
                arrowprops=dict(arrowstyle="->"))
    a0.set_xlabel("physical error probability $p$ (depolarising, per data qubit)"); a0.set_ylabel("logical error probability")
    a0.set_title("(a) Code capacity: perfect syndrome extraction"); a0.legend(fontsize=7); a0.grid(alpha=.3, which="both")
    a1 = ax[1]
    for (k, v), c in zip(circ.items(), ["C1", "C2", "C4"]):
        a1.errorbar(cps, [x[0] for x in v], yerr=[x[1] for x in v], fmt="o-", color=c, capsize=2, label=k)
    a1.errorbar(cps, [x[0] for x in bare], yerr=[x[1] for x in bare], fmt="ks--", capsize=2, label="bare physical qubit")
    a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlabel("physical error probability $p$ (every gate and measurement)")
    a1.set_ylabel("logical failure probability"); a1.set_title("(b) Circuit-level noise: non-fault-tolerant circuits")
    a1.legend(fontsize=7); a1.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "main_result.png"), dpi=200); plt.close(fig)

    # weight figure
    fc = S.failing_counts(); fig, ax2 = plt.subplots(figsize=(4.6, 3))
    ax2.bar(list(fc), [fc[w][0] / fc[w][1] for w in fc], color="C0"); ax2.set_xlabel("weight of Pauli error")
    ax2.set_ylabel("fraction causing a logical error"); ax2.set_title("Where the code breaks")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "failure_by_weight.png"), dpi=200); plt.close(fig)

    out["runtime_s"] = round(time.time() - t0)
    json.dump(out, open(os.path.join(ROOT, "results.json"), "w"), indent=1, default=str)
    print("finished in", out["runtime_s"], "s")


if __name__ == "__main__":
    main()
