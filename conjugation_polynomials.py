from sage.all import var, PolynomialRing, GF, factor, Matrix, vector, QQ

x = var("x")
P = PolynomialRing(GF(2), x)
Q = PolynomialRing(QQ, x)

for i in range(1, 2**7):
    init = [int(a) for a in bin(i)[2:].zfill(8)]
    num = vector(init) * vector([x**i for i in range(8)])
    den = 1 - x**8
    t = (P(num) / P(den)).denominator()
    n = t.degree()
    M = Matrix([[i**j for j in range(n)] for i in range(n)])
    s = M.solve_right(vector(init[:n]))
    res = s * vector([x**i for i in range(n)])
    if Q(res).degree() <= 4:
        print(init, factor(res))
