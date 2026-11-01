import unittest

from simple_chinese_remainder_theorem import crt, CRTError


class TestCRT(unittest.TestCase):
    def test_single_congruence(self):
        x, M = crt([3], [7])
        self.assertEqual((x, M), (3, 7))

    def test_classic_sunzi(self):
        # x == 2 (mod 3), x == 3 (mod 5), x == 2 (mod 7)
        # The smallest positive solution is 23.
        x, M = crt([2, 3, 2], [3, 5, 7])
        self.assertEqual((x, M), (23, 105))

    def test_two_coprime_moduli(self):
        x, M = crt([1, 2], [3, 5])
        self.assertEqual((x, M), (7, 15))

    def test_residues_larger_than_moduli_are_reduced(self):
        x, M = crt([10, 12], [3, 5])
        # 10 % 3 == 1, 12 % 5 == 2, which is the same system as (1, 2).
        self.assertEqual((x, M), (7, 15))

    def test_zero_residues(self):
        x, M = crt([0, 0, 0], [2, 3, 5])
        self.assertEqual((x, M), (0, 30))

    def test_modulus_one_is_allowed(self):
        # 1 is coprime to everything, and a single congruence modulo 1
        # forces no constraint, so the answer should match the other one.
        x, M = crt([5, 0], [7, 1])
        self.assertEqual((x, M), (5, 7))

    def test_solution_is_canonical(self):
        x, M = crt([2, 3, 2], [3, 5, 7])
        self.assertTrue(0 <= x < M)

    def test_coprimality_of_many_moduli(self):
        moduli = [2, 3, 5, 7, 11, 13]
        residues = [1, 2, 3, 4, 5, 6]
        x, M = crt(residues, moduli)
        self.assertEqual(M, 2 * 3 * 5 * 7 * 11 * 13)
        for r, m in zip(residues, moduli):
            self.assertEqual(x % m, r)

    def test_large_moduli(self):
        p = 10**9 + 7
        q = 10**9 + 9
        x, M = crt([123456789, 987654321], [p, q])
        self.assertEqual(M, p * q)
        self.assertEqual(x % p, 123456789)
        self.assertEqual(x % q, 987654321)

    def test_empty_input_raises(self):
        with self.assertRaises(CRTError):
            crt([], [])

    def test_mismatched_lengths_raise(self):
        with self.assertRaises(CRTError):
            crt([1, 2], [3])

    def test_non_positive_modulus_raises(self):
        with self.assertRaises(CRTError):
            crt([1], [0])
        with self.assertRaises(CRTError):
            crt([1], [-5])

    def test_non_coprime_moduli_raise(self):
        with self.assertRaises(CRTError):
            crt([1, 2], [4, 6])

    def test_non_coprime_in_later_pair_raises(self):
        # 2 and 4 are coprime to 9, but 2 and 4 share a factor.
        with self.assertRaises(CRTError):
            crt([1, 2, 3], [2, 4, 9])

    def test_negative_residues_are_normalised(self):
        x, M = crt([-1, -2], [5, 7])
        # -1 % 5 == 4 and -2 % 7 == 5.
        self.assertEqual((x, M), (19, 35))


if __name__ == "__main__":
    unittest.main()
