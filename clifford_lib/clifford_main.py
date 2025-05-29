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


if __name__ == "__main__":
    p = 7
    q = 0
    r = 0
    Cl = get_clifford(p, q, r)
    gens = [val for _, val in Cl.basis(1).items()]
