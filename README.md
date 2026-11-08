# Coincidence site lattices of the cubic lattice, up to symmetry

The number of coincidence site lattices (CSLs) of the simple cubic lattice Z^3 of each index, when
lattices that are related by a rotation of the cube are counted as the same.

A CSL is what two copies of a crystal lattice have in common after one copy is rotated. For special
rotations a fraction 1/Σ of the lattice points coincide, and Σ is called the index. CSLs are the
standard way of describing grain boundaries in cubic crystals, which is where the familiar Σ3, Σ5,
Σ7, ... boundaries come from.

The total number of CSLs of index 2n-1 is OEIS [A031359](https://oeis.org/A031359). Many of those
lattices are just rotated copies of each other. This repository counts the genuinely different ones.

## The sequence

a(n) for index Σ = 2n-1, n = 1, 2, 3, ...

```
1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 2, 3, 4, 3, 3, 5, 3, ...
```

For example, there are 14 CSLs of index 13, and they fall into two classes under the rotations of
the cube, of sizes 6 and 8, so a(7) = 2. The 48 CSLs of index 33 fall into three classes, of sizes
12, 12 and 24, so a(17) = 3.

Counting under all 48 symmetries of the cube instead of the 24 rotations gives the same numbers,
since every lattice is mapped to itself by x -> -x.

## Files

- `csl_orbits.py` builds every coincidence rotation from a primitive integer quaternion, computes
  each CSL exactly as the dual of Z^3 + RZ^3, and groups the lattices into orbits under the
  symmetries of the cube. `python csl_orbits.py 99` covers every odd index up to 99.
- `csl_orbits_snf.py` is a second program written separately. It uses the same rotations but gets
  each CSL from the Smith normal form of the rotation matrix. `python csl_orbits_snf.py 41`.
- `csl.gp` is a PARI/GP program for the closed formula (Zeiner 2005, Theorem 2.2), plus a
  function that writes `terms.txt`.
- `terms.txt` lists n and a(n) for n = 1 to 10000, from `csl.gp`.

Both Python programs use the standard library only.

## Checks

- Both programs find exactly A031359(n) lattices of each index, so no lattice is missed or
  counted twice.
- The class counts from `csl_orbits.py` agree with the formula for every odd index up to 99
  (n = 1 to 50), and `csl_orbits_snf.py` agrees for every odd index up to 41 (n = 1 to 21).

## References

- H. Grimmer, The generating function for coincidence site lattices in the cubic system,
  Acta Cryst. A40 (1984), 108-112.
- P. Zeiner, Symmetries of coincidence site lattices of cubic lattices, Z. Kristallogr. 220
  (2005), 915-925. [arXiv:math/0605525](https://arxiv.org/abs/math/0605525)
- M. Baake, Solution of the coincidence problem in dimensions d <= 4, in R. V. Moody (ed.),
  The Mathematics of Long-Range Aperiodic Order, Kluwer, 1997, 9-44.
  [arXiv:math/0605222](https://arxiv.org/abs/math/0605222)

## License

MIT
