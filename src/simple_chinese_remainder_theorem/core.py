"""Chinese Remainder Theorem for pairwise-coprime moduli."""

from __future__ import annotations

from typing import List, Sequence, Tuple


class CRTError(ValueError):
    """Raised when the inputs cannot be solved by the CRT."""


def _extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Return (g, x, y) such that a*x + b*y = g = gcd(a, b).

    Iterative implementation so that large moduli do not hit Python's
    recursion limit.  The sign of g is always non-negative.
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def _mod_inverse(a: int, m: int) -> int:
    """Return x in [0, m) with a*x == 1 (mod m).

    Raises CRTError if a and m are not coprime.
    """
    if m <= 0:
        raise CRTError("modulus must be positive")
    g, x, _ = _extended_gcd(a % m, m)
    if g != 1:
        raise CRTError(f"{a} has no inverse modulo {m}: gcd is {g}")
    return x % m


def crt(residues: Sequence[int], moduli: Sequence[int]) -> Tuple[int, int]:
    """Solve a system of congruences with pairwise-coprime moduli.

    Given ``residues`` and ``moduli`` of equal length, find the unique
    integer ``x`` modulo ``M = product(moduli)`` such that
    ``x % moduli[i] == residues[i] % moduli[i]`` for every ``i``.

    Returns ``(x, M)`` with ``0 <= x < M``.

    The implementation is the standard incremental CRT: fold the
    congruences left to right, combining the running solution with the
    next one using the extended Euclidean algorithm.  This keeps the
    intermediate values bounded by the running product, so it never
    needs to compute the full product before starting.

    Raises CRTError if:
      * the two sequences have different lengths;
      * the input is empty (the empty system has no unique solution);
      * any modulus is not strictly positive;
      * any two moduli are not coprime (the system may then have zero or
        many solutions, and this function only handles the unique case).
    """
    if len(residues) != len(moduli):
        raise CRTError("residues and moduli must have the same length")
    if not moduli:
        raise CRTError("cannot solve an empty system")

    # Validate all moduli up front so the error message names the bad one
    # rather than failing midway through the fold.
    for i, m in enumerate(moduli):
        if m <= 0:
            raise CRTError(f"modulus at index {i} must be positive, got {m}")

    # Check pairwise coprimality explicitly.  The incremental algorithm
    # would also catch it via _mod_inverse, but the message would be
    # confusing (it would mention a non-invertible residue).
    for i in range(len(moduli)):
        for j in range(i + 1, len(moduli)):
            g, _, _ = _extended_gcd(moduli[i], moduli[j])
            if g != 1:
                raise CRTError(
                    f"moduli at indices {i} and {j} are not coprime: "
                    f"gcd({moduli[i]}, {moduli[j]}) = {g}"
                )

    x = residues[0] % moduli[0]
    M = moduli[0]
    for i in range(1, len(moduli)):
        r = residues[i] % moduli[i]
        m = moduli[i]
        # Solve x + M*t == r (mod m), i.e. M*t == r - x (mod m).
        diff = (r - x) % m
        inv = _mod_inverse(M, m)
        t = (diff * inv) % m
        x = x + M * t
        M *= m
        x %= M

    return x, M
