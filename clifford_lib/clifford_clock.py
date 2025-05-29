import itertools as it
import functools as ft
import random

from sage.all import QQ, CliffordAlgebra
from sage.all import I as i
from sage.all import (
    Matrix,
    PolynomialRing,
    QuadraticForm,
    QuaternionAlgebra,
    block_diagonal_matrix,
    block_matrix,
    identity_matrix,
    prod,
    vector,
    zero_matrix,
)
from tqdm import tqdm

from clifford_main import (
    get_clifford,
    get_clifford_random_element_from_gens,
    get_full_clifford_from_gens,
    get_full_random_from_gens,
)

R = QQ

# Qx = PolynomialRing(QQ, "x")
# x = Qx.gen()
# RR = Qx.quotient(x**2 - 1)
# C = Qx.quotient(x**2 + 1)
# RR = MatrixSpace(R, 2)
# C = MatrixSpace(R, 2)

# H = QuaternionAlgebra(-1, -1)
# Hx = PolynomialRing(H, "x")
# x = Hx.gen()
# HH = Hx.quotient(x**2 + 1)
# HH = MatrixSpace(H, 2)


def to_mat(x, basis_to_mat: dict, inv_basis: dict):
    res = 0
    for el, idx in inv_basis.items():
        res += basis_to_mat[el] * x[idx]
    return res


def from_mat(y, basis_to_mat: dict, basis: list):
    n_unknowns = len(basis_to_mat)
    assert (n_unknowns & (n_unknowns - 1)) == 0
    n_equations = prod(basis_to_mat[basis[0]].dimensions())

    def unroll_matrix(M):
        return sum([list(row) for row in M], [])

    res_matrix = Matrix(n_equations, n_unknowns)
    for i, e in enumerate(basis):
        column = Matrix(unroll_matrix(basis_to_mat[e])).T
        res_matrix.set_block(0, i, column)

    res_vector = vector(unroll_matrix(y))

    # print(res_matrix)
    # print(res_vector)

    assert res_matrix.rank() == n_unknowns
    ans = res_matrix.solve_right(res_vector)
    return sum(x * y for x, y in zip(ans, basis))


def base_case(p, q, gens: list):
    gens_to_mat = {}
    n = p + q

    Hgen1 = Matrix(R, [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
    Hgen2 = Matrix(R, [[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]])
    Hgen3 = Matrix(R, [[0, 0, 0, -1], [0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]) # - k = i * j

    if (p, q) == (0, 0):
        pass

    elif (p, q) == (1, 0):
        gens_to_mat[gens[0]] = Matrix(R, [[0, 1], [1, 0]])

    elif (p, q) == (0, 1):
        gens_to_mat[gens[0]] = Matrix(R, [[0, 1], [-1, 0]])

    elif (p, q) == (0, 2):
        gens_to_mat[gens[0]] = Hgen1
        gens_to_mat[gens[1]] = Hgen2
        # gens[inv_basis[basis[3]]] = Hgen3

    # u = ue + u1e1 + u2e2 + u12e12 + (u123 - u23e1 + u13e2 - u3e12) * (e123 = omega, omega^2 = 1, omega commutes)
    # u = ue + u1e1 + u2e2 + u12e12 + u3e3 + u13e13 + u23e23 + u123e123
    elif (p, q) == (0, 3):
        gens_to_mat[gens[0]] = block_diagonal_matrix(Hgen1, Hgen1)
        gens_to_mat[gens[1]] = block_diagonal_matrix(Hgen2, Hgen2)
        # gens[inv_basis[basis[4]]] = block_diagonal_matrix(Hgen3, Hgen3)

        Zero = zero_matrix(R, 4)
        gens_to_mat[gens[2]] = block_matrix(2, 2, [Zero, -Hgen3, -Hgen3, Zero])

        # Id = identity_matrix(R, 4)
        # gens[inv_basis[basis[5]]] = block_matrix(2, 2, [Zero, Hgen2, Hgen2, Zero])
        # gens[inv_basis[basis[6]]] = block_matrix(2, 2, [Zero, -Hgen1, -Hgen1, Zero])
        # gens[inv_basis[basis[7]]] = block_matrix(2, 2, [Zero, Id, Id, Zero])
    else:
        raise ValueError(f"Wrong base case: {p}, {q}")

    return gens_to_mat


def internal_get_isomorphism(p, q, gensp: list, gensq: list):
    # print(p, q)
    n = p + q
    sig = p - q

    if (p, q) in [(0, 0), (0, 1), (1, 0), (0, 2), (0, 3)]:
        return base_case(p, q, gensp + gensq)

    elif sig >= 2:
        gensp_new = [gensp[0]] + [eq * gensp[0] for eq in gensq]
        gensq_new = [ep * gensp[0] for ep in gensp[1:]]

        gens_to_mat = internal_get_isomorphism(q + 1, p - 1, gensp_new, gensq_new)

        new_gens_to_mat = {}

        # 1st unchanged
        new_gens_to_mat[gensp[0]] = gens_to_mat[gensp_new[0]]

        gtm0 = gens_to_mat[gensp_new[0]]
        # gtm0 = gens_to_mat[gensp_new[0]]**-1
        # print("get rid of this after")
        # assert gtm0 == gens_to_mat[gensp_new[0]]

        # 2..p from q = p - 1 els
        for i, e in enumerate(gensq_new):
            new_gens_to_mat[gensp[i + 1]] = gens_to_mat[e] * gtm0

        # 1...q from p - 1 = q + 1 - 1 = q
        for i, e in enumerate(gensp_new[1:]):
            new_gens_to_mat[gensq[i]] = gens_to_mat[e] * gtm0

        assert len(new_gens_to_mat) == p + q
        return new_gens_to_mat

    elif sig == -n:
        w = prod(gensq[:4])
        gensp_new = [e * w for e in gensq[:4]] + gensp.copy()
        gensq_new = gensq.copy()[4:]
        gens_to_mat = internal_get_isomorphism(p + 4, q - 4, gensp_new, gensq_new)

        new_gens_to_mat = {}

        w1 = prod(gensp_new[:4])
        W1 = prod(gens_to_mat[e] for e in gensp_new[:4])

        assert gensp_new[0] * w1 == gensq[0]
        assert gensp_new[1] * w1 == gensq[1]
        assert gensp_new[2] * w1 == gensq[2]
        assert gensp_new[3] * w1 == gensq[3]

        new_gens_to_mat[gensq[0]] = gens_to_mat[gensp_new[0]] * W1
        new_gens_to_mat[gensq[1]] = gens_to_mat[gensp_new[1]] * W1
        new_gens_to_mat[gensq[2]] = gens_to_mat[gensp_new[2]] * W1
        new_gens_to_mat[gensq[3]] = gens_to_mat[gensp_new[3]] * W1

        for i, e in enumerate(gensq_new):
            new_gens_to_mat[gensq[i + 4]] = gens_to_mat[e]

        for i, e in enumerate(gensp_new[4:]):
            new_gens_to_mat[gensp[i]] = gens_to_mat[e]

        assert len(new_gens_to_mat) == p + q
        return new_gens_to_mat

    else:
        gensp_new = gensp[:-1]
        gensq_new = gensq[:-1]
        gens_to_mat = internal_get_isomorphism(p - 1, q - 1, gensp_new, gensq_new)

        new_gens_to_mat = {}

        for i, e in enumerate(gensp_new):
            new_gens_to_mat[gensp[i]] = block_diagonal_matrix(
                gens_to_mat[e], -gens_to_mat[e]
            )

        for i, e in enumerate(gensq_new):
            new_gens_to_mat[gensq[i]] = block_diagonal_matrix(
                gens_to_mat[e], -gens_to_mat[e]
            )

        if len(new_gens_to_mat) == 0:
            dim = 1
        else:
            dim = list(gens_to_mat.values())[0].dimensions()[0]

        Zero = zero_matrix(R, dim)
        Id = identity_matrix(R, dim)
        new_gens_to_mat[gensp[-1]] = block_matrix(2, 2, [Zero, Id, Id, Zero])
        new_gens_to_mat[gensq[-1]] = block_matrix(2, 2, [Zero, -Id, Id, Zero])
        assert len(new_gens_to_mat) == p + q
        return new_gens_to_mat


def recover_basis(Cl, gens, gens_to_mat):
    n = len(gens)

    basis = [Cl(1)] + gens.copy()
    basis_to_mat = gens_to_mat.copy()

    dim = basis_to_mat[basis[1]].dimensions()[0]
    basis_to_mat[Cl(1)] = identity_matrix(R, dim)

    for i in range(2, n + 1):
        for e in it.combinations(gens, r=i):
            ng = prod(e)
            basis.append(ng)
            basis_to_mat[ng] = prod([gens_to_mat[e_i] for e_i in e])

    return basis, basis_to_mat


def get_isomorphism(p, q):
    n = p + q
    Cl = get_clifford(p, q)
    gens = [v for _, v in Cl.basis(1).items()]
    gensp = gens[:p]
    gensq = gens[p:]

    gens_to_mat = internal_get_isomorphism(p, q, gensp, gensq)

    basis, basis_to_mat = recover_basis(Cl, gens, gens_to_mat)
    inv_basis = {j: i for k in range(n + 1) for i, j in Cl.basis(k).items()}

    iso = ft.partial(to_mat, basis_to_mat=basis_to_mat, inv_basis=inv_basis)
    iso_inv = ft.partial(from_mat, basis_to_mat=basis_to_mat, basis=basis)
    return Cl, basis, iso, iso_inv


def check_isomorphism(p, q):
    Cl, basis, iso, iso_inv = get_isomorphism(p, q)
    gens = basis[1: p + q + 1]

    u = Cl(get_full_random_from_gens(gens))
    v = Cl(get_full_random_from_gens(gens))

    u_m = iso(u)
    v_m = iso(v)

    lhs = u * v
    rhs = iso_inv(u_m * v_m)
    assert lhs == rhs

    lhs = iso(u * v)
    rhs = u_m * v_m
    assert lhs == rhs

def check():
    for i, j in it.product(range(5), repeat=2):
        if i + j == 0:
            continue
        check_isomorphism(i, j)

    check_isomorphism(4, 5)

if __name__ == "__main__":
    check()

#p, q = 4, 5
#Cl, basis, basis_to_mat, inv_basis, gens = get_isomorphism(p, q)
#u = Cl(get_full_random_from_gens(gens))
#u_m = to_mat(u, basis_to_mat, inv_basis)
#print(u_m)
