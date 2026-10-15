"""Run the illustrative, non-quantum Steane syndrome arithmetic demo."""
from .steane_syndrome import check_single


def main() -> None:
    passed = sum(check_single(position, pauli)
                 for position in range(7) for pauli in ("X", "Y", "Z"))
    print(f"{passed}/21 algebraic single-qubit Pauli checks passed")
    print("NOTE: This is NOT a Qiskit circuit or hardware experiment.")
    if passed != 21:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
