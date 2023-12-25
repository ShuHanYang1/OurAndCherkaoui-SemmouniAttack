import time


# n is the size of primes, alpha is the size of e, MSB_size ist the bit-size of the unknown MSBs.
mm=8
print('m = ', mm)
n = 100
alpha = 400
#betla 表示的是 p q 相同的部分 位置的部分是1/2- betla
betla=0.125
detla=0.72
nbetla = int(n * betla * 2)
rnbetla = int(n - nbetla)
print("nbetla=",nbetla)
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


#生成p。q低位相同
def keyGen(n,alpha, nbetla,rnbetla):

    while (1):

        if nbetla == 0:
            S = 0
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



    print("p=",p)
    print("bin(p)=",bin(p))
    print("q=",q)
    print("bin(q)=", bin(q))
    print("e=", e)
    print("d=", d)



    return N, p, q, e,d,S,k


N, p, q, e, d,S,k= keyGen(n, alpha, nbetla,rnbetla)
start_time = time.time()

PR.< u,x, y > = PolynomialRing(ZZ)

XX = int(3 * N ^ (alpha / (2 * n) + detla - 2))
YY = int(2 * N ^ (1  - 2*betla))

A = -(N - 1) ^ 2
pol = 1 + A * x + x * y
Q = PR.quotient(x*y +1 - u) # u = xy + 1
polZ = Q(pol).lift()
UU = XX*YY + 1
print(polZ)
re=(p-q)^2

# print((x*polZ).monomials())


gg = []
for kk in range(mm + 1):
    for ii in range(mm - kk + 1):
        xshift = x^ ii * e ^ (mm - kk) * polZ(u, x, y) ^ kk
        #print("ii,kk:",ii,kk,( x^ ii * e ^ (mm - kk) * polZ(u, x, y) ^ kk ).monomials())
        gg.append(xshift)
gg.sort()

# x-shifts list of monomials
monomials = []
for polynomial in gg:
    for monomial in polynomial.monomials():
        if monomial not in monomials:
            monomials.append(monomial)
monomials.sort()

tt =(1-detla+2*betla)*mm/(1-2*betla)
# print("tao=",int((1-detla+2*betla)/(1-2*betla)))
print("tt=",tt)


# y-shifts (selected by Herrman and May)
for jj in range(1, int(tt) + 1):
    #for kk in range(0, mm + 1):
    for kk in range(floor(mm/ tt)*jj , mm + 1):
        #print("jj",jj,"kk",kk)
        yshift = y ^ jj * polZ(u, x, y) ^ kk * e ^ (mm - kk)
        yshift = Q(yshift).lift()
        gg.append(yshift)  # substitution

# y-shifts list of monomials
for jj in range(1, int(tt) + 1):
    #for kk in range( mm + 1):
    for kk in range(floor(mm / tt) * jj, mm + 1):
        monomials.append(u ^ kk * y ^ jj)

nnm = len(monomials)
# for ii in range(nnm):
#     print(monomials[ii],end=' ')
# print()

# construct lattice B
nn = len(monomials)
BB = Matrix(ZZ, nn)
for ii in range(nn):
    BB[ii, 0] = gg[ii](0, 0, 0)
    for jj in range(1, ii + 1):
        if monomials[jj] in gg[ii].monomials():
            BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](UU, XX, YY)

matrix_overview(BB,e ^ mm)


det=BB.det()
bound = e ^ (mm * nn)
if abs(det) >= bound:
    print("We do not have det < bound. Solutions might not be found.")
    print("Try with highers m and t.")
    diff = (log(det) - log(bound)) / log(2)
    print(log(det*1.0)/log(2))
    #print("size det(L) - size e^(m*n) = ", floor(diff))
else:
    print("det(L) < e^(m*n) (good! If a solution exists < N^delta, it will be found)")


BB = BB.LLL()


matrix_overview(BB, e ^ mm)
print("Orgin mm=",mm," detla=",detla," n=",n)
for pol1_idx in range(nn - 1):
    for pol2_idx in range(pol1_idx + 1, nn):
        # for i and j, create the two polynomials
        PR.< w, z > = PolynomialRing(ZZ)
        pol1 = pol2 = 0
        for jj in range(nn):
            pol1 += monomials[jj](w * z + 1, w, z) * BB[pol1_idx, jj] / monomials[jj](UU, XX, YY)
            pol2 += monomials[jj](w * z + 1, w, z) * BB[pol2_idx, jj] / monomials[jj](UU, XX, YY)

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
print("x==k:",solx==-k,solx)
print("y==-(p-q)^2:",soly==re,re,soly)
print("=== %sLLL seconds ===" % (time.time() - start_time))








