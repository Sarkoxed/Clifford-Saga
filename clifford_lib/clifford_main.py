import itertools as it
import random

from sage.all import QQ, CliffordAlgebra, QuadraticForm, binomial, prod
from tqdm import tqdm


def get_clifford(p, q, r=0):
    n = p + q + r

    Mq = []
    for i in range(p):
        Mq.extend([1] + [0] * (n - i - 1))
    for i in range(p, p + q):
        Mq.extend([-1] + [0] * (n - i - 1))
    for i in range(p + q, n):
        Mq.extend([0] * (n - i))

    Qf = QuadraticForm(QQ, n, Mq)
    Cl = CliffordAlgebra(Qf, "e")
    return Cl


def get_clifford_random_element_from_gens(gens, S=20):
    N = len(gens)
    start = random.randint(-10, 10)
    while start in QQ:
        for i in range(random.randint(1, S)):
            k = random.randint(0, N)
            start += random.randint(-10, 10) * prod(random.choices(gens, k=k))
    return start


def get_full_clifford_from_gens(gens, Cl):
    N = len(gens)
    start = 1
    for i in range(1, N + 1):
        for j in it.combinations(gens, r=i):
            start += prod(j)
    return start


def get_full_random_from_gens(gens):
    N = len(gens)
    start = random.randint(-20, 20)
    for i in range(1, N + 1):
        for j in it.combinations(gens, r=i):
            start += random.randint(-20, 20) * prod(j)
    return start


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


if __name__ == "__main__":
    p = 7
    q = 0
    r = 0
    Cl = get_clifford(p, q, r)
    gens = [val for _, val in Cl.basis(1).items()]

# for _, e in Cl.basis():
#        e = prod(k)
#        print(e, e**2)
# exit()

id = lambda x: 0

## reverse
reverse = lambda x: binomial(x, 2)
# for i in tqdm(range(100)):
#    u = get_clifford_random_element(gens)
#    v = get_clifford_random_element(gens)
#    lhs = abstract_conj(reverse, u * v, N)
#    rhs = abstract_conj(reverse, v, N) * abstract_conj(reverse, u, N)
#    assert lhs == rhs
#
## oddness
oddness = lambda x: x
# for i in tqdm(range(100)):
#    u = get_clifford_random_element(gens)
#    v = get_clifford_random_element(gens)
#    lhs = abstract_conj(oddness, u * v, N)
#    rhs = abstract_conj(oddness, u, N) * abstract_conj(oddness, v, N)
#    assert lhs == rhs
#
## under
under = lambda x: x != 0
# for i in tqdm(range(100)):
#    u = get_clifford_random_element(gens)
#    v = get_clifford_random_element(gens)
#    lhs = abstract_conj(under, u * v, N) * u
#    rhs = u * abstract_conj(under, v * u, N)
#    assert lhs == rhs

tetrate = lambda x: binomial(x, 4)

clifford = lambda x: binomial(x, 2) + x

# print(get_full_clifford(gens))

# f = get_full_clifford(gens)
# print(f)
# for k in get_mods(f, N, 4):
#    print(k)
# hs = get_mods(f, N, r)
# h0, h1, h2, h3 = hs
# print(hs)


# r = 8
# u = get_full_random(gens)
# v = get_full_random(gens)
# us = get_mods(u, N, r)
# vs = get_mods(v, N, r)
#
# pairs = [(i, i) for i in range(r)] + list(it.combinations(range(r), r=2))
# for i, j in pairs:
#    x, y = us[i], vs[j]
#    x = abstract_conj(clifford, x, N)
#    t = x * y - y * x
#    # print(i, j, x * y - y * x)
#    print(i, j, get_rank(t))
