"""Two independent evaluations of the signed moment matrices at both embeddings."""
from mpmath import mp

mp.dps = 165
N = 265

def lucas_values(P, Q, n):
    values = [0, 1]
    for _ in range(2, n+1):
        values.append(P*values[-1]-Q*values[-2])
    return values

def formal_moments(param, Q, size, X):
    rho, w = Q*param, Q*param**2
    roots = [w**(-a) for a in range(1,size+1)]
    pweights = []
    ds = [mp.mpf(1)]
    S = 0
    for a in range(1,size+1):
        wa = w**a
        S += param**(3*a)/(1-wa)**3
        partial = mp.fprod(1-w**(b-a) for b in range(1,size+1) if a!=b)
        pweights.append(2*rho**(-a)*(X-S)/partial)
        ds = [v-wa*u for v,u in zip(ds+[0],[0]+ds)]
    values=[mp.fsum(c*x**j for c,x in zip(pweights,roots)) for j in range(size)]
    for j in range(size-1):
        alpha=rho*w**j
        nu=j*(j+1)*alpha/(1-alpha)
        values.append((nu-mp.fsum(d*values[j+h] for h,d in enumerate(ds[:-1])))/ds[-1])
    return values

def direct_real(param, Q, size):
    # Apply (s*d/ds)^2 + s*d/ds directly to s^k / D(s).
    rho,w=Q*param,Q*param**2
    roots=[w**a for a in range(1,size+1)]
    values=[mp.mpf(0) for _ in range(2*size-1)]
    for n in range(1,N+1):
        s = w**n
        factors=[s*r for r in roots]
        first=mp.fsum(f/(1-f) for f in factors)
        second=mp.fsum(f/(1-f)**2 for f in factors)
        weight=rho**n/mp.fprod(1-f for f in factors)
        for k in range(len(values)):
            values[k]+=weight*((k+first)**2+(k+first)+second)
            weight*=s
    return values

def contour_coefficients(param, Q, size):
    rho,w=Q*param,Q*param**2
    # At each r_a the conjugated simple-pole moment equals Gamma(r_a).
    # Expand Gamma and 1/F on the unit circle to extract the contour coefficient.
    roots=[w**a for a in range(1,size+1)]
    h=[mp.mpf(1)]+[mp.mpf(0)]*(N+size)
    for r in roots:
        for j in range(1,len(h)):
            h[j]+=r*h[j-1]
    nu=[]; nu_conj=[]
    for j in range(N+1):
        alpha=rho*w**j
        nu.append(j*(j+1)*alpha/(1-alpha))
        nu_conj.append(Q*j*(j+1)/(alpha-Q))
    values=[]
    for k in range(2*size-1):
        proper=mp.fsum(nu[j]*h[k+j-size+1] for j in range(max(1,size-1-k),N+1))
        quotient=mp.fsum(h[l]*nu_conj[k-size-l] for l in range(k-size+1))
        scale=(-1)**size*w**(size*(size+1)//2)
        values.append((proper+quotient)*scale)
    return values

def relative(x,y):
    return abs(x-y)/max(abs(x),abs(y),mp.mpf('1e-100'))

def hankel(values,K):
    return mp.det(mp.matrix([[values[i+j] for j in range(K)] for i in range(K)]))

def one(P,Q):
    u_seq=lucas_values(P,Q,N)
    D=P*P-4*Q
    root=mp.sqrt(D)
    param=-2/(P+root)
    w=Q*param**2
    A=mp.fsum((-1)**(n-1)/mp.mpf(u_seq[n])**3 for n in range(1,N+1))
    Xi=mp.fsum(param**(3*n)/(1-w**n)**3 for n in range(1,N+1))
    assert relative(A,-D*root*Xi)<mp.mpf('1e-145')
    max_error=mp.mpf(0)
    for K in (2,3,4,5):
        direct=direct_real(param,Q,K)
        formal=formal_moments(param,Q,K,Xi)
        conjugate=formal_moments(Q/param,Q,K,-Xi)
        analytic=contour_coefficients(param,Q,K)
        for a,b in [(formal,direct),(conjugate,analytic)]:
            err=max(relative(x,y) for x,y in zip(a,b))
            max_error=max(max_error,err,relative(hankel(a,K),hankel(b,K)))
            assert hankel(a,K)!=0
    assert max_error<mp.mpf('1e-88'),max_error
    negatives=lucas_values(-P,Q,40)
    assert all(negatives[n]==(-1)**(n-1)*u_seq[n] for n in range(1,len(negatives)))
    traces=[2,P]
    for k in range(2,9):
        traces.append(P*traces[-1]-Q*traces[-2])
    for k in range(1,9):
        new=lucas_values(traces[k],Q**k,24)
        for n in range(1,25):
            assert u_seq[k*n]==u_seq[k]*new[n]
    print(f"(P,Q)=({P},{Q:+d}), A={mp.nstr(A,19)}, K=2..5, maximum relative discrepancy {mp.nstr(max_error,3)}",flush=True)

if __name__=='__main__':
    for P,Q in [(1,-1),(2,-1),(3,1)]:
        one(P,Q)
