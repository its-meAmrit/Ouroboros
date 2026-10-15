"""Steane [[7,1,3]] syndrome *algebra*, not a quantum circuit simulation.

H columns are all seven distinct nonzero three-bit strings. For a Pauli
error, its X and Z binary components give the Z-check and X-check syndromes.
"""

H_COLUMNS = (1, 2, 3, 4, 5, 6, 7)


def syndrome(bitmask: int) -> int:
    """Return the three-bit Hamming syndrome of a seven-bit error mask."""
    if not 0 <= bitmask < 128:
        raise ValueError("bitmask must be between 0 and 127")
    result = 0
    for i, col in enumerate(H_COLUMNS):
        if bitmask & (1 << i):
            result ^= col
    return result


def pauli_components(position: int, pauli: str) -> tuple[int, int]:
    """Return X, Z component masks; Y has both components."""
    if position not in range(7) or pauli not in {"X", "Y", "Z"}:
        raise ValueError("position must be 0..6 and pauli X, Y, or Z")
    bit = 1 << position
    return (bit if pauli in {"X", "Y"} else 0,
            bit if pauli in {"Z", "Y"} else 0)


def decode_single(syndrome_x: int, syndrome_z: int) -> tuple[int, int]:
    """Decode assumed single-qubit X and Z components from both syndromes.

    `syndrome_x` is detected by Z checks; `syndrome_z` by X checks.
    Distinct nonzero syndromes may correspond to multi-qubit errors: this
    lookup only promises single-error correction, not fault tolerance.
    """
    lookup = {col: 1 << i for i, col in enumerate(H_COLUMNS)}
    if syndrome_x not in range(8) or syndrome_z not in range(8):
        raise ValueError("syndromes must be 0..7")
    return lookup.get(syndrome_x, 0), lookup.get(syndrome_z, 0)


def check_single(position: int, pauli: str) -> bool:
    x_mask, z_mask = pauli_components(position, pauli)
    corrected_x, corrected_z = decode_single(syndrome(x_mask), syndrome(z_mask))
    return (x_mask ^ corrected_x, z_mask ^ corrected_z) == (0, 0)
