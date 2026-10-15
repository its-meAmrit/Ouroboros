"""Check toy algebra only; not encoded-state correctness."""
import unittest
from src.steane_syndrome import syndrome, pauli_components, check_single, decode_single


class SteaneSyndromeAlgebraTests(unittest.TestCase):
    def test_each_single_pauli(self):
        for i in range(7):
            for p in ("X", "Y", "Z"):
                with self.subTest(qubit=i, pauli=p):
                    self.assertTrue(check_single(i, p))

    def test_no_error(self):
        self.assertEqual(syndrome(0), 0)
        self.assertEqual(decode_single(0, 0), (0, 0))

    def test_each_column_distinct(self):
        self.assertEqual({syndrome(1 << i) for i in range(7)}, set(range(1, 8)))

    def test_y_has_both_components(self):
        self.assertEqual(pauli_components(2, "Y"), (4, 4))

    def test_invalid_position(self):
        with self.assertRaises(ValueError):
            pauli_components(7, "X")


if __name__ == "__main__":
    unittest.main()
