"""Exact small determinants in Q(T)[X], and the cost of their cyclotomic divisor."""
from fractions import Fraction
from math import prod
from sympy import QQ
from sympy.polys.fields import field

field_t, t = field("t", QQ)
ring = field_t.ring
T = ring.gens[0]

def cluster_cost(m, singular_first):
    stop = m if singular_first else m - 1
    return sum(max(0, m + 2 - 2*j) for j in range(stop))

def clearing_exponents(K, Q):
    shift = Fraction(int(Q == -1), 2)
    result = {}
    for d in range(1, 4*K + 1):
        order = (Fraction(2, d) + shift).denominator
        cost = 0
        if order <= K:
            for r in range(1, order+1):
                size = 1 + (K-r)//order
                cost += cluster_cost(size, r == order)
        for j in range(1, K-1):
            alpha = Fraction(2*j+1, d) + (j+1)*shift
            if alpha.denominator == 1:
                cost += K-1-j
                break
        if cost:
            result[d] = cost
    return result

def cyclotomic_polynomials(limit):
    polys = {}
    for d in range(1, limit+1):
        lower = prod((polys[h] for h in range(1, d) if d%h == 0), start=ring.one)
        polys[d] = (T**d-1).exquo(lower)
    return polys

def entries(K, Q):
    w, rho = Q*t**2, Q*t
    nodes = [w**(-a) for a in range(1, K+1)]
    denom = [field_t.one]
    simple_values = []
    tail_removed = field_t.zero
    for a in range(1, K+1):
        root_factor = w**a
        denom = [v-root_factor*u for v,u in zip(denom+[0],[0]+denom)]
        tail_removed += t**(3*a)/(1-w**a)**3
        residue = 1 / prod((1-w**(b-a) for b in range(1,K+1) if b != a), start=field_t.one)
        slope = residue*2*rho**(-a)
        simple_values.append((-tail_removed*slope,slope))
    out = []
    for k in range(K):
        out.append(tuple(sum((simple_values[j][s]*nodes[j]**k for j in range(K)),start=field_t.zero) for s in (0,1)))
    for j in range(K-1):
        alpha = rho*w**j
        moment = j*(j+1)*alpha/(1-alpha)
        out.append(tuple(((moment if s == 0 else 0)-sum(denom[h]*out[j+h][s] for h in range(K)))/denom[-1] for s in (0,1)))
    return out

def determinant(matrix):
    n = len(matrix)
    sign = 1
    previous = ring.one
    for k in range(n-1):
        pivot_row = next((i for i in range(k,n) if matrix[i][k]),None)
        if pivot_row is None:
            return ring.zero
        if pivot_row != k:
            matrix[k],matrix[pivot_row] = matrix[pivot_row],matrix[k]
            sign = -sign
        pivot = matrix[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                matrix[i][j]=(pivot*matrix[i][j]-matrix[i][k]*matrix[k][j]).exquo(previous)
        previous = pivot
    return sign*matrix[-1][-1]

def lagrange_coefficients(value, K):
    coeffs=[1]
    div=1
    for other in range(K+1):
        if other == value:
            continue
        coeffs=[a-other*b for a,b in zip([0]+coeffs,coeffs+[0])]
        div*=value-other
    return [QQ(c,div) for c in coeffs]

def check(K,Q,cyclos):
    exponents = clearing_exponents(K,Q)
    C = prod((cyclos[d]**e for d,e in exponents.items()),start=ring.one)
    pairs = entries(K,Q)
    common = ring.one
    for pair in pairs:
        for f in pair:
            common=common.lcm(f.denom)
    pairs = [tuple(f.numer*common.exquo(f.denom) for f in pair) for pair in pairs]
    delta = [ring.zero for _ in range(K+1)]
    for x in range(K+1):
        seq=[c+x*m for c,m in pairs]
        value=determinant([[seq[i+j] for j in range(K)] for i in range(K)])
        for i,weight in enumerate(lagrange_coefficients(x,K)):
            delta[i]+=weight*value
    assert delta[-1], "degree for the tested small size"
    actual = ring.one
    for coef in delta:
        raw = field_t(coef)/field_t(common**K)
        actual = actual.lcm(raw.denom)
        cleared = field_t(C)*raw
        assert len(cleared.denom) == 1
        scale = next(iter(cleared.denom.values()))
        assert all((v/scale).denominator == 1 for v in cleared.numer.values())
    actual = actual.monic()
    while actual and not actual.get((0,),0):
        actual=actual.exquo(T)
    assert C.rem(actual) == 0
    print(f"Q={Q:+d} K={K} degree in X={K}, true cyclotomic degree={actual.degree()}, allowed degree={C.degree()}",flush=True)

def cubic_constants():
    def weight(h):
        return prod((1-Fraction(1,p) for p in (2,3,5,7) if h%p == 0),start=Fraction(1))
    beta=sum(weight(h) for h in range(1,211))/210
    assert beta == prod(1-Fraction(1,p*p) for p in (2,3,5,7))
    # Integral from 0 to 1 of 8*x*(1-x) + max(0,1-4*x)^2/2.
    w_integral=8*(Fraction(1,2)-Fraction(1,3))+Fraction(1,24)
    cost=w_integral*beta+Fraction(1,3)
    margin=1+Fraction(2,9)-cost
    assert (beta,cost,margin)==(Fraction(768,1225),Fraction(4393,3675),Fraction(296,11025))
    print(f"sieve mean {beta}, integral {w_integral}, cost {cost}, margin {margin}",flush=True)

if __name__=="__main__":
    cubic_constants()
    cyclo=cyclotomic_polynomials(20)
    for k in (2,3,4,5):
        for q in (-1,1):
            check(k,q,cyclo)
