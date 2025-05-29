from clifford_lib.clifford_clock import get_isomorphism
from clifford_lib.clifford_main import get_full_random_from_gens

p, q = 4, 5
Cl, basis, iso, iso_inv = get_isomorphism(p, q)
gens = basis[1 : p + q + 1]

u = Cl(get_full_random_from_gens(gens))
u_m = iso(u)
print(u_m)
