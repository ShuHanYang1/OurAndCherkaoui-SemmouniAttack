import time
# n is the size of primes, alpha is the size of e, MSB_size ist the bit-size of the unknown MSBs.
mm=4

print('m = ', mm)
n = 100
alpha = 400
#betla 表示的是 p q 相同的部分 位置的部分是1/2- betla
betla=0.24
detla=1.4
nbetla = int(n * betla * 2)
rnbetla = int(n - nbetla)
print(nbetla)
print(rnbetla)
print('alpha= ',alpha/(2*n) ,' & betla= ',betla,' & detla= ',detla)


def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        num = 0
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            if BB[ii, jj] != 0:
                num = jj
            a += '0' if BB[ii, jj] == 0 else 'X'
            if BB.dimensions()[0] < 60:
                a += ' '
        if abs(BB[ii, num]) >= bound:
            a += '~'
        print(a)


def matrix(m,tt,w,alpha,betla,delta,N):
    detx=1
    dety=1
    det=1
    for k in range(0,m+1):
        for l in range(0,m-k):
            detx*=pow(N,(alpha/(2*n))*m+((alpha/(2*n))+delta-2)*l+(delta-1-4*betla)*k)
    for i in range(2, floor(tt * m) + 1):
        for j in range(floor(i / tt), m + 1):
            dety*=pow(N,(alpha/(2*n))*m+(1/2-2*betla)*i+(delta-1-4*betla)*j)

    A = pow(((m + 1) * floor(tt * m-1)),(w/2)) * pow((1 + pow(m, 3 * m)), pow(w, 2))
    det=A*detx*dety
    det1=detx*dety
    return det,det1


#生成p。q低位相同
def keyGen(n,alpha, nbetla,rnbetla):

    while (1):
        if nbetla==0:
            S=0
        else:
            S = ZZ.random_element(2 ^ (nbetla - 1), 2 ^ (nbetla))
        pL = next_prime(ZZ.random_element(2 ^ (rnbetla)))
        qL = next_prime(ZZ.random_element(2 ^ (rnbetla)))


        p =S * 2 ^ (rnbetla) + pL
        q = S * 2 ^ (rnbetla) + qL

        if is_prime(p) == False or is_prime(q) == False:
            continue

        ndetla = int(2*n * detla)
        d = next_prime(ZZ.random_element(2 ^ (ndetla - 1), 2 ^ (ndetla)))
        N = p * q

        phi = (p ^ 2 - 1) * (q ^ 2 - 1)
        if (N.nbits() != 2 * n or gcd(d, (p ^ 2 - 1) * (q ^ 2 - 1)) != 1):
            continue

        e = d.inverse_mod(phi)
        k=(e*d-1)/phi
        if (e.nbits() == alpha):  # N=pq should be 2n bits
            break
        #print(k)



    print("p=",p)
    print("bin(p)=",bin(p))
    print("q=",q)
    print("bin(q)=", bin(q))
    print("e=", e)
    print("d=", d)



    return N, p, q, e,d,S,k


N, p, q, e, d,S,k= keyGen(n, alpha, nbetla,rnbetla)

def getV(N,e,S,p,q):
    if S==0:
        return 0
    M=floor(N/(2^(2*n-nbetla)))
    M=floor(sqrt(M*2^nbetla))


    v=M*(2^nbetla)+(floor((N - M ^ 2 * 2 ^ (2 * (rnbetla))) / (M * 2 ^ (rnbetla))) - floor(
        (N - M ^ 2 * 2 ^ (2 * (rnbetla))) / (M * 2 ^ (rnbetla))) % 2 ^ (n - 2 * nbetla + 1)) / 2 ^ (n - 2 * nbetla + 1)

    sum=p+q
    vT=(sum-sum%2^(n-2*nbetla+1))/2^(n-2*nbetla+1)


    return v

start_time = time.time()

v=getV(N,e,S,p,q)


PR.<u, x, y > = PolynomialRing(ZZ)
XX = int(3 * N ^ (alpha / (2 * n) + detla - 2))
YY = int(2 * N ^ (1 / 2 - 2 * betla))
UU = XX*YY^2 - 1
A = (N + 1) ^ 2 - (v * 2 ^ (n - 2 * nbetla + 1)) ^ 2
a = v * 2 ^ (n - 2 * nbetla + 2)
polZ = 1 + A * x - a * x * y - x * y ^ 2
Q = PR.quotient(x*y^2 - 1 - u) # u = xy^2 - 1
pol = Q(polZ).lift()

res=(p+q)%2^(n - 2 * nbetla + 1)



# x-shifts
gg = []
monomials = []
detch=1
for kk in range(mm + 1):
    for ii in range(mm - kk + 1):
        for jj in range(2):
            detch=detch*XX^(ii+kk)*YY^(2*kk+jj)*e^(mm-kk)
            xshift = x ^ ii * y ^ jj * pol(u,x, y) ^ kk * e ^ (mm - kk)
            xshift = Q(xshift).lift()
            for monomial in xshift.monomials():
                if monomial not in monomials:
                    monomials.append(monomial)
            gg.append(xshift)


tt = (2 * (1 + 4 * betla - detla))/(1 - 4 * betla)
print("tt=",tt)

# y-shifts (selected by Herrman and May)
ymonomials = []
ynum=0
for jj in range(2,floor(tt*mm) + 1):
    for kk in range(floor(jj/tt),mm+1):
        ynum=ynum+1
        yshift = y ^ jj * pol(u,x, y) ^ kk * e ^ (mm - kk)
        yshift = Q(yshift).lift()
        for monomial in yshift.monomials():
            if monomial not in ymonomials:
                ymonomials.append(monomial)
        gg.append(yshift)


detch=detch*pow(((mm + 1) * floor(tt * mm-1)),(ynum/2)) * pow((1 + pow(mm, 3 * mm)), pow(ynum, 2))
#ymonomials.sort()
for monomial in ymonomials:
    if monomial not in monomials:
        monomials.append(monomial)

# construct lattice B
nnm = len(monomials)

nng = len(gg)
BB = Matrix(ZZ, nng, nnm)


for ii in range(nng):
    BB[ii, 0] = gg[ii](0, 0,0)
    for jj in range(1, nnm):
        if monomials[jj] in gg[ii].monomials():
            BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](UU,XX, YY)

matrix_overview(BB,e ^ mm)

CC = BB * (BB.transpose())
det = int(sqrt(CC.det()))




bound = e ^ (mm * nng)




if abs(det) >= bound:
    print("We do not have det < bound. Solutions might not be found.")
    print("Try with highers m and t.")
    diff = (log(det) - log(bound)) / log(2)
    #print(diff)
else:
    print("det(L) < e^(m*n) (good! If a solution exists < N^delta, it will be found)")



BB = BB.LLL()


matrix_overview(BB, e ^ mm)

found_polynomials = False

for pol1_idx in range(nng - 1):
    for pol2_idx in range(pol1_idx + 1, nng):
        # for i and j, create the two polynomials
        PR.< w, z > = PolynomialRing(ZZ)
        pol1 = pol2 = 0
        for jj in range(nnm):
            pol1 += monomials[jj](w * z^2 - 1, w, z) * BB[pol1_idx, jj] / monomials[jj](UU,XX, YY)
            pol2 += monomials[jj](w * z^2 - 1, w, z) * BB[pol2_idx, jj] / monomials[jj](UU,XX, YY)

        # resultant
        PR.< q > = PolynomialRing(ZZ)
        rr = pol1.resultant(pol2)

        # are these good polynomials?
        if rr.is_zero() or rr.monomials() == [1]:
            continue
        else:
            print("found them, using vectors", pol1_idx, "and", pol2_idx)
            found_polynomials = True
            break
    if found_polynomials:
        break

if not found_polynomials:
    print("no independant vectors could be found. This should very rarely happen...")


rr = rr(q, q)

# solutions
soly = rr.roots()

if len(soly) == 0:
    print("Your prediction (delta) is too small")


soly = soly[0][0]
ss = pol1(q, soly)
solx = ss.roots()[0][0]
print('x = ', solx==k,solx,k)
print('y = ',res==soly,soly,res)

print("=== %all seconds ===" % (time.time() - start_time))












