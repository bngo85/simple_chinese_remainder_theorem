# simple_chinese_remainder_theorem

Solve a system of congruences with pairwise-coprime moduli and return the unique solution modulo their product, using the extended Euclidean algorithm.

```python
from simple_chinese_remainder_theorem import crt, CRTError

# x == 2 (mod 3), x == 3 (mod 5), x == 2 (mod 7)
x, M = crt([2, 3, 2], [3, 5, 7])
print(x, M)  # 23 105

# Residues are reduced modulo their modulus, so these are equivalent:
print(crt([10, 12], [3, 5]))  # (7, 15)

try:
    crt([1, 2], [4, 6])  # 4 and 6 are not coprime
except CRTError as exc:
    print(exc)  # moduli at indices 0 and 1 are not coprime: gcd(4, 6) = 2
```

The two exported names are `crt` and `CRTError`.

## Why this exists

Python's standard library has no built-in Chinese Remainder Theorem solver, and the CRT is the natural way to combine results that were computed modulo several pairwise-coprime numbers — for example, reconstructing a value from residues, or checking a computation against several moduli. The trade-off made here is strictness: the function only accepts pairwise-coprime moduli, and it raises `CRTError` otherwise. That keeps the contract simple and the result genuinely unique. A more permissive solver would have to decide what to do when the system is inconsistent or has multiple solutions, and those decisions are better made by the caller.

## The awkward edge

Residues are normalised with `%` before use, so negative residues and residues larger than their modulus are accepted. The returned solution is always in `[0, M)`. If you pass a modulus that is not strictly positive, or two moduli that share a factor, `crt` raises `CRTError` rather than guessing. An empty input also raises `CRTError`, because the empty system has no unique solution.
