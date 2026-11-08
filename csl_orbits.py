"""Coincidence site lattices (CSLs) of the cubic lattice Z^3, counted up to the symmetries of the cube.

Every coincidence rotation of Z^3 has the form R(q) = M(q)/|q|^2 for a primitive integer quaternion q.
For each such rotation the CSL is Z^3 intersected with R Z^3. Since Z^3 is its own dual and R is orthogonal,
this intersection is the dual of the sum Z^3 + R Z^3, which is computed exactly in integer arithmetic.
Each CSL is stored by its Hermite normal form, and the CSLs of each index are grouped into orbits under the
48 signed permutation matrices.

For a primitive quaternion, |q|^2 is never divisible by 8, so the index |q|^2 / 2^k is at least |q|^2 / 4.
Searching all quaternions with |q|^2 <= 4 * SMAX therefore finds every CSL with index up to SMAX.

Usage: python csl_orbits.py [SMAX]      (default 99; only odd indices occur)
"""
import itertools
import sys
from fractions import Fraction
from math import gcd, isqrt


def hnf(vectors):
    """Lower triangular Hermite normal form of the lattice spanned by integer vectors in Z^3 (full rank)."""
    vecs = [list(v) for v in vectors if any(v)]
    basis = []
    for r in range(3):
        while True:
            nonzero = [v for v in vecs if v[r] != 0]
            if len(nonzero) <= 1:
                break
            nonzero.sort(key=lambda v: abs(v[r]))
            p = nonzero[0]
            for v in nonzero[1:]:
                q = v[r] // p[r]
                for k in range(3):
                    v[k] -= q * p[k]
            vecs = [v for v in vecs if any(v)]
        pivot = next(v for v in vecs if v[r] != 0)
        if pivot[r] < 0:
            pivot[:] = [-x for x in pivot]
        vecs = [v for v in vecs if v is not pivot]
        basis.append(pivot)
    for c in range(3):
        for r in range(c + 1, 3):
            q = basis[c][r] // basis[r][r]
            for k in range(3):
                basis[c][k] -= q * basis[r][k]
    return tuple(tuple(v) for v in basis)


def rotation(q):
    a, b, c, d = q
    m = a * a + b * b + c * c + d * d
    M = [[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
         [2*(b*c + a*d), a*a - b*b + c*c - d*d, 2*(c*d - a*b)],
         [2*(b*d - a*c), 2*(c*d + a*b), a*a - b*b - c*c + d*d]]
    return [[Fraction(x, m) for x in row] for row in M]


def csl(R):
    """Return (index, HNF of the CSL) for the rational rotation R."""
    d = 1
    for row in R:
        for x in row:
            d = d * x.denominator // gcd(d, x.denominator)
    # columns of d * (Z^3 + R Z^3)
    gens = [[d if i == j else 0 for i in range(3)] for j in range(3)]
    gens += [[int(R[i][j] * d) for i in range(3)] for j in range(3)]
    B = hnf(gens)
    Bm = [[Fraction(B[j][i]) for j in range(3)] for i in range(3)]
    det = (Bm[0][0] * (Bm[1][1] * Bm[2][2] - Bm[1][2] * Bm[2][1])
           - Bm[0][1] * (Bm[1][0] * Bm[2][2] - Bm[1][2] * Bm[2][0])
           + Bm[0][2] * (Bm[1][0] * Bm[2][1] - Bm[1][1] * Bm[2][0]))
    cof = [[Bm[(i+1) % 3][(j+1) % 3] * Bm[(i+2) % 3][(j+2) % 3]
            - Bm[(i+1) % 3][(j+2) % 3] * Bm[(i+2) % 3][(j+1) % 3] for j in range(3)] for i in range(3)]
    # dual basis of (Z^3 + R Z^3) is d * B^(-T); B^(-T) = cof / det
    dual = [[cof[i][j] / det * d for j in range(3)] for i in range(3)]
    assert all(x.denominator == 1 for row in dual for x in row)
    cols = [[int(dual[i][j]) for i in range(3)] for j in range(3)]
    L = hnf(cols)
    return abs(L[0][0] * L[1][1] * L[2][2]), L


def symmetries():
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            out.append((perm, signs))
    return out


def transform(g, L):
    perm, signs = g
    return hnf([[signs[i] * v[perm[i]] for i in range(3)] for v in L])


def psi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            r = r // p * (p + 1)
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        r = r // m * (m + 1)
    return r


def main():
    smax = int(sys.argv[1]) if len(sys.argv) > 1 else 99
    bound = isqrt(4 * smax) + 1
    G = symmetries()
    found = {}
    for q in itertools.product(range(-bound, bound + 1), repeat=4):
        if sum(x * x for x in q) > 4 * smax or q == (0, 0, 0, 0):
            continue
        if gcd(gcd(q[0], q[1]), gcd(q[2], q[3])) != 1:
            continue
        s, L = csl(rotation(q))
        if s <= smax:
            found.setdefault(s, set()).add(L)
    print("index  CSLs  psi  orbits  orbit sizes")
    counts = []
    for s in range(1, smax + 1, 2):
        lattices = found.get(s, set())
        orbits = {}
        for L in lattices:
            orbits.setdefault(min(transform(g, L) for g in G), []).append(L)
        sizes = sorted(len(v) for v in orbits.values())
        ok = "" if len(lattices) == psi(s) else "  (total differs from psi!)"
        print(f"{s:5}  {len(lattices):4}  {psi(s):3}  {len(orbits):6}  {sizes}{ok}")
        counts.append(len(orbits))
    print()
    print("a(n) for n = 1..%d:" % len(counts))
    print(", ".join(map(str, counts)))


if __name__ == "__main__":
    main()
