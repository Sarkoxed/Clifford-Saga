from sage.all import QQ, CliffordAlgebra, QuadraticForm, binomial, prod

from .clifford_main import get_clifford_random_element_from_gens


def get_mods(t, N, M):
    if t in QQ:
        return [t] + [0] * (M - 1)
    mods = []
    for i in range(M):
        tmp = 0
        for j in range(i, N + 1, M):
            tmp += t.homogeneous_component(j)
        mods.append(tmp)
    return mods


def get_rank(t):
    s = [a.degree() for a in t.monomials()]
    res = list(set(s))
    if len(res) == 0:
        res = [0]
    return res


def abstract_conj(f, t, n):
    if t in QQ:
        return t
    res = 0
    for i in range(n + 1):
        res += (-1) ** f(i) * t.homogeneous_component(i)
    return res


id = lambda x: 0
reverse = lambda x: binomial(x, 2)
oddness = lambda x: x
under = lambda x: x != 0
tetrate = lambda x: binomial(x, 4)
clifford = lambda x: binomial(x, 2) + x

if __name__ == "__main__":
    for i in tqdm(range(100)):
        u = get_clifford_random_element_from_gens(gens)
        v = get_clifford_random_element_from_gens(gens)
        lhs = abstract_conj(reverse, u * v, N)
        rhs = abstract_conj(reverse, v, N) * abstract_conj(reverse, u, N)
        assert lhs == rhs

    for i in tqdm(range(100)):
        u = get_clifford_random_element_from_gens(gens)
        v = get_clifford_random_element_from_gens(gens)
        lhs = abstract_conj(oddness, u * v, N)
        rhs = abstract_conj(oddness, u, N) * abstract_conj(oddness, v, N)
        assert lhs == rhs

    for i in tqdm(range(100)):
        u = get_clifford_random_element_from_gens(gens)
        v = get_clifford_random_element_from_gens(gens)
        lhs = abstract_conj(under, u * v, N) * u
        rhs = u * abstract_conj(under, v * u, N)
        assert lhs == rhs
