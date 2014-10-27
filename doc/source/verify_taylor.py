from gmpy import mpq, lcm, denom, numer, fac

def atan_coefficients(NN, bits):
    ps = []
    qs = []
    temp = []
    Q = 1
    for k in range(2*NN+50):
        p = 1
        q = 2*k+1
        if lcm(Q, q) < 2**bits:
            temp.append(mpq(p,q))
            Q = lcm(Q, q)
        else:
            for a in temp:
                ps.append(int(a * Q))
                qs.append(int(Q))
            Q = q
            temp = [mpq(p,q)]
    return ps[:NN], qs[:NN]

def exp_coefficients(M, bits):
    N = 2*M+50
    Qs = [fac(k) for k in range(N)]
    prevstop = 0
    for k in range(N):
        if Qs[k] >= 2**bits-1:
            q = Qs[k-1]
            for i in range(k, N): Qs[i] //= q
            for i in range(prevstop, k): Qs[i] = q
            prevstop = k
    Ps = Qs[:]
    fact = 1
    for k in range(1, N):
        assert Qs[k] < 2**bits-1
        if Qs[k] == Qs[k-1]:
            fact *= k
        else:
            fact = k
        Ps[k] //= fact
    return map(int, Ps)[:N], map(int, Qs)[:N]

class FixedPointBound(object):

    def __init__(self, bits, mid, rad):
        self.bits = bits
        self.mid = mpq(mid)
        self.rad = mpq(rad)  # rad is in ulp

    def add(self, other):
        if isinstance(other, FixedPointBound):
            mid = self.mid + other.mid
            rad = self.rad + other.rad
        else:
            assert other == int(other) and other >= 0
            mid = self.mid + int(other)
            rad = self.rad
        return FixedPointBound(self.bits, mid, rad)

    def mul(self, other):
        if isinstance(other, FixedPointBound):
            MAX_ULP = mpq(1, 2**self.bits)
            mid = self.mid * other.mid
            rad = 0
            rad += self.rad * other.mid  # ulp
            rad += self.mid * other.rad  # ulp
            rad += self.rad * other.rad * MAX_ULP  # ulp
            rad += 1   # ulp rounding
        else:
            assert other == int(other) and other >= 0
            mid = self.mid * int(other)
            rad = self.rad * int(other)
        return FixedPointBound(self.bits, mid, rad)

    def div(self, other):
        assert other == int(other) and other >= 0
        mid = self.mid / mpq(other)
        rad = self.rad / mpq(other) + 1
        return FixedPointBound(self.bits, mid, rad)

    def addmul(self, other, c):
        assert c == int(c) and c >= 0
        c = abs(int(c))
        mid = self.mid + other.mid * c
        rad = self.rad + other.rad * c
        return FixedPointBound(self.bits, mid, rad)

    def check_overflow_0(self):
        # check that self fits 0 integral limbs
        MAX_ULP = mpq(1, 2**self.bits)
        assert self.mid + self.rad * MAX_ULP < 1 - MAX_ULP

    def check_overflow_1(self):
        # check that self fits 1 integral limb
        MAX_ULP = mpq(1, 2**self.bits)
        assert self.mid + self.rad * MAX_ULP < 2**self.bits - MAX_ULP

    def check_le_int(self, c):
        # check that |self| <= c
        MAX_ULP = mpq(1, 2**self.bits)
        assert self.mid + self.rad * MAX_ULP <= c

def verify_atan(N, PS, QS, bits):
    X = FixedPointBound(bits, mpq(1,16), 0)
    S = FixedPointBound(bits, 0, 0)
    m = 2
    while m * m < N:
        m += 2
    T = [None] * (m+1)
    T[1] = X.mul(X)
    T[2] = T[1].mul(T[1])
    for k in range(4, m + 1, 2):
        T[k-1] = T[k//2].mul(T[k//2-1])
        T[k] = T[k//2].mul(T[k//2])
    for k in range(N-1, -1, -1):
        c, d, e = PS[k], QS[k], QS[k+1]
        if d != e and k < N-1:
            # if alternating, adding e must give a nonnegative number
            S.check_le_int(e)
            # adding e must not overflow
            S.add(e).check_overflow_1()
            S = S.mul(d).div(e)
            # if alternating, adding d must not overflow
            S.add(d).check_overflow_1()
        if k % m == 0:
            # if alternating, adding c must give a nonnegative number
            S.check_le_int(c)
            S = S.add(c)
            S.check_overflow_1()
            if k != 0:
                S = S.mul(T[m])
                S.check_overflow_1()
        else:
            S = S.addmul(T[k % m], c)
            S.check_overflow_1()
    S = S.div(mpq(QS[0]))
    S = S.mul(X)
    S.check_overflow_0()
    print N, float(S.mid), float(S.rad)
    assert S.rad <= 2

def verify_exp(N, PS, QS, bits):
    X = FixedPointBound(bits, mpq(1,16), 0)
    S = FixedPointBound(bits, 0, 0)
    m = 2
    while m * m < N:
        m += 2
    T = [None] * (m+1)
    T[1] = X
    T[2] = T[1].mul(T[1])
    for k in range(4, m + 1, 2):
        T[k-1] = T[k//2].mul(T[k//2-1])
        T[k] = T[k//2].mul(T[k//2])
    for k in range(N-1, -1, -1):
        c, d, e = PS[k], QS[k], QS[k+1]
        if d != e and k < N-1:
            # if alternating, adding e must give a nonnegative number
            S.check_le_int(e)
            # adding e must not overflow
            S.add(e).check_overflow_1()
            S = S.div(e)
            # if alternating, adding 1 must not overflow
            S.add(1).check_overflow_1()
        if k % m == 0:
            # if alternating, adding c must give a nonnegative number
            S.check_le_int(c)
            S = S.add(c)
            S.check_overflow_1()
            if k != 0:
                S = S.mul(T[m])
                S.check_overflow_1()
        else:
            S = S.addmul(T[k % m], c)
            S.check_overflow_1()
    S = S.div(mpq(QS[0]))
    S.check_overflow_1()
    print N, float(S.mid), float(S.rad)
    assert S.rad <= 2

def verify_sin_cos(N, PS, QS, bits):
    X = FixedPointBound(bits, mpq(1,16), 0)
    m = 2
    while m * m < N:
        m += 2
    T = [None] * (m+1)
    T[1] = X.mul(X)
    T[2] = T[1].mul(T[1])
    for k in range(4, m + 1, 2):
        T[k-1] = T[k//2].mul(T[k//2-1])
        T[k] = T[k//2].mul(T[k//2])
    for cosorsin in range(2):
        S = FixedPointBound(bits, 0, 0)
        for k in range(N-1, -1, -1):
            c, d, e = PS[2*k+cosorsin], QS[2*k+cosorsin], QS[2*k+cosorsin+2]
            if d != e and k < N-1:
                # if alternating, adding e must give a nonnegative number
                S.check_le_int(e)
                # adding e must not overflow
                S.add(e).check_overflow_1()
                S = S.div(e)
                # if alternating, adding 1 must not overflow
                S.add(1).check_overflow_1()
            if k % m == 0:
                # if alternating, adding c must give a nonnegative number
                S.check_le_int(c)
                S = S.add(c)
                S.check_overflow_1()
                if k != 0:
                    S = S.mul(T[m])
                    S.check_overflow_1()
            else:
                S = S.addmul(T[k % m], c)
                S.check_overflow_1()
        if cosorsin == 0:
            S = S.div(mpq(QS[0]))
            S.check_overflow_1()
            # note: top limb must actually be 0 or 1;
            # but this follows by S.rad <= 2
            print N, float(S.mid), float(S.rad)
            assert S.rad <= 2
        else:
            S = S.div(mpq(QS[0]))
            S.check_overflow_1()
            S = S.mul(X)
            S.check_overflow_0()
            print N, float(S.mid), float(S.rad)
            assert S.rad <= 2

for bits in [32, 64]:
    PS, QS = exp_coefficients(300, bits)
    for N in range(300):
        verify_sin_cos(N, PS, QS, bits)

for bits in [32, 64]:
    PS, QS = exp_coefficients(300, bits)
    for N in range(300):
        verify_exp(N, PS, QS, bits)

for bits in [32, 64]:
    PS, QS = atan_coefficients(300, bits)
    for N in range(300):
        verify_atan(N, PS, QS, bits)

