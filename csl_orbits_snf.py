"""Second, independent count of the cubic coincidence site lattices up to the symmetries of the cube.

Same quaternion parametrisation of the coincidence rotations as csl_orbits.py, but the lattice is built
differently. For q = (a,b,c,d) with N = a^2 + b^2 + c^2 + d^2 and integer matrix M = N * R(q),

    CSL = { x in Z^3 : M x = 0 (mod N) }.

If M = U D V is the Smith normal form (D diagonal, U and V unimodular), then M x = 0 (mod N) exactly when
the i-th coordinate of V x is divisible by e_i = N / gcd(D_ii, N). So the CSL is V^(-1) diag(e) Z^3, with
index e_1 e_2 e_3. Lattices are identified by Hermite normal form and grouped into orbits under the 48
signed permutation matrices. Every rotation is checked to be orthogonal with determinant +1.

Usage: python csl_orbits_snf.py [SMAX]      (default 41)
"""
import sys
from itertools import permutations, product
from math import gcd, isqrt


def smith_3(Min):
    """Diagonal of the Smith normal form of a 3x3 integer matrix, and the right factor V (M = U D V)."""
    M = [row[:] for row in Min]
    n = 3
    V = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def col_op(j1, j2, q):          # column j1 -= q * column j2, mirrored on the rows of V
        for i in range(n):
            M[i][j1] -= q * M[i][j2]
        for k in range(n):
            V[j2][k] += q * V[j1][k]

    def col_swap(j1, j2):
        for i in range(n):
            M[i][j1], M[i][j2] = M[i][j2], M[i][j1]
        V[j1], V[j2] = V[j2], V[j1]

    def row_op(i1, i2, q):
        for j in range(n):
            M[i1][j] -= q * M[i2][j]

    def row_swap(i1, i2):
        M[i1], M[i2] = M[i2], M[i1]

    for t in range(n):
        while True:
            piv = None
            for i in range(t, n):
                for j in range(t, n):
                    if M[i][j] != 0:
                        piv = (i, j)
                        break
                if piv:
                    break
            if piv is None:
                break
            if piv[0] != t:
                row_swap(t, piv[0])
            if piv[1] != t:
                col_swap(t, piv[1])
            done = True
            for i in range(t + 1, n):
                if M[i][t] != 0:
                    row_op(i, t, M[i][t] // M[t][t])
                    if M[i][t] != 0:
                        row_swap(t, i)
                        done = False
            for j in range(t + 1, n):
                if M[t][j] != 0:
                    col_op(j, t, M[t][j] // M[t][t])
                    if M[t][j] != 0:
                        col_swap(t, j)
                        done = False
            if done:
                break
    return [M[i][i] for i in range(n)], V


def inv_unimodular(V):
    a, b, c = V[0]
    d, e, f = V[1]
    g, h, i = V[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    adj = [[(e * i - f * h), -(b * i - c * h), (b * f - c * e)],
           [-(d * i - f * g), (a * i - c * g), -(a * f - c * d)],
           [(d * h - e * g), -(a * h - b * g), (a * e - b * d)]]
    assert det in (1, -1), det
    return [[adj[r][c2] // det for c2 in range(3)] for r in range(3)]


def hnf_columns(cols):
    """Lower triangular Hermite normal form of three independent integer columns."""
    A = [[cols[j][i] for j in range(3)] for i in range(3)]
    for i in range(3):
        while True:
            nz = [j for j in range(i, 3) if A[i][j] != 0]
            if not nz:
                break
            j0 = min(nz, key=lambda j: abs(A[i][j]))
            if j0 != i:
                for r in range(3):
                    A[r][i], A[r][j0] = A[r][j0], A[r][i]
            if A[i][i] == 0:
                break
            cleared = True
            for j in range(i + 1, 3):
                if A[i][j] != 0:
                    q = A[i][j] // A[i][i]
                    for r in range(3):
                        A[r][j] -= q * A[r][i]
                    if A[i][j] != 0:
                        cleared = False
            if cleared:
                break
        if A[i][i] < 0:
            for r in range(3):
                A[r][i] = -A[r][i]
    for i in range(3):
        for j in range(i):
            if A[i][i] != 0:
                q = A[i][j] // A[i][i]
                for r in range(3):
                    A[r][j] -= q * A[r][i]
    return tuple(tuple(row) for row in A)


def rot_matrix(a, b, c, d):
    return [[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]]


def check_rotation(M, N):
    for i in range(3):
        for j in range(3):
            if sum(M[k][i] * M[k][j] for k in range(3)) != (N * N if i == j else 0):
                return False
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
           - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    return det == N ** 3


PERMS = []
for perm in permutations(range(3)):
    for signs in product((1, -1), repeat=3):
        m = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        for i in range(3):
            m[i][perm[i]] = signs[i]
        PERMS.append(m)


def apply_mat(P, cols):
    return [tuple(sum(P[i][k] * col[k] for k in range(3)) for i in range(3)) for col in cols]


def psi(m):
    out, x, p = 1, m, 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out *= (p + 1) * p ** (e - 1)
        p += 1
    if x > 1:
        out *= x + 1
    return out


def run(sigma_max):
    qbound = isqrt(4 * sigma_max) + 1
    by_index = {}
    checked = 0
    for a in range(0, qbound + 1):
        for b in range(-qbound, qbound + 1):
            for c in range(-qbound, qbound + 1):
                for d in range(-qbound, qbound + 1):
                    if a == 0 and (b, c, d) <= (0, 0, 0):
                        continue
                    N = a * a + b * b + c * c + d * d
                    if N == 0 or N > 4 * sigma_max:
                        continue
                    if gcd(gcd(a, b), gcd(c, d)) != 1:
                        continue
                    M = rot_matrix(a, b, c, d)
                    if not check_rotation(M, N):
                        raise SystemExit("rotation check failed")
                    checked += 1
                    D, V = smith_3(M)
                    e = [N // gcd(abs(dd) if dd else N, N) for dd in D]
                    idx = e[0] * e[1] * e[2]
                    if idx > sigma_max or idx % 2 == 0:
                        continue
                    Vi = inv_unimodular(V)
                    cols = [tuple(Vi[r][j] * e[j] for r in range(3)) for j in range(3)]
                    by_index.setdefault(idx, set()).add(hnf_columns(cols))

    print("quaternions checked: %d" % checked)
    print("index  CSLs  psi  orbits")
    orbit_counts = []
    for s in range(1, sigma_max + 1, 2):
        lats = by_index.get(s, set())
        seen, orbits = set(), 0
        for L in sorted(lats):
            if L in seen:
                continue
            orbits += 1
            cols = [tuple(L[r][j] for r in range(3)) for j in range(3)]
            for P in PERMS:
                seen.add(hnf_columns(apply_mat(P, cols)))
        orbit_counts.append(orbits)
        print("%5d  %4d  %3d  %6d" % (s, len(lats), psi(s), orbits))
        if len(lats) != psi(s):
            raise SystemExit("index %d: found %d lattices, expected psi = %d" % (s, len(lats), psi(s)))
    print()
    print("a(n) for n = 1..%d:" % len(orbit_counts))
    print(", ".join(map(str, orbit_counts)))


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 41)
