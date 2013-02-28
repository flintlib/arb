**gamma.h** -- support for the gamma function
===============================================================================

This module implements various algorithms for evaluating the
gamma function and related functions. The functions provided here are mainly
intended for internal use, though they may be useful to call directly in some
applications where the default algorithm choices are suboptimal.
Most applications should use the standard, user-friendly top-level functions:

* :func:`fmprb_gamma`
* :func:`fmprb_rgamma`
* :func:`fmprb_lgamma`
* :func:`fmpcb_gamma`
* :func:`fmpcb_rgamma`
* :func:`fmpcb_lgamma`


Evaluation using Taylor series
--------------------------------------------------------------------------------

.. function :: long gamma_taylor_bound_mag(long n)

    Letting `1/\Gamma(x) = \sum_{n=1}^{\infty} c_n x^n`, returns an integer
    `p` such that `|c_n| \le 2^{-p}`.
    Uses the bound `\log |c_n| \le b(n) = ((\pi-1)n + (3-5n) \log(n/6))/6`
    (to be published).

.. function :: void gamma_taylor_bound_ratio(fmpr_t r, long n)

    Sets `r` to a bound for `B(n+1) / B(n)` where `B(n) = \exp b(n)`.
    We use `r(n) = 3n^{-5/6}`.

.. function :: void gamma_taylor_bound_remainder(fmpr_t err, const fmpr_t z, long n)

    Given a nonnegative `z` and `n \ge 1`, computes a bound for
    `\sum_{k=n}^{\infty} |c_{k+1}| z^k`. This is given by
    `B(n+1) z^n / (1 - r(n+1) z)`, provided that the denominator is
    positive (otherwise the bound is infinity).

.. function :: void gamma_taylor_fmprb(fmprb_t y, const fmprb_t x, long prec)

    Sets `y = \Gamma(x)`, computed by translating `x` to the interval
    `[-1/2,1/2]` and then evaluating the Taylor series for the
    reciprocal gamma function.

    Assumes that the absolute value of the midpoint of `x` is not too large
    in order for the translation to be reasonably efficient (too large
    values will  result in extreme slowness and eventually an exception).

Evaluation using the Stirling series
--------------------------------------------------------------------------------

.. function :: void gamma_stirling_choose_param(int * reflect, long * r, long * n, double x, double y, double beta, int allow_reflection, long prec)

    Uses double precision arithmetic to compute parameters `r`, `n` such that
    the remainder in the Stirling series with `z = x+yi`
    approximately is bounded by `2^{-\mathrm{prec}}`.

    The parameter `n` is the truncation point in the asymptotic
    Stirling series. If `|z|` is too small for the Stirling series
    to give sufficient accuracy directly, we first translate to `z + r`
    using the formula `\Gamma(z) = \Gamma(z+r) / 
    (z (z+1) (z+2) \cdots (z+r-1))`.

    If *allow_reflection* is nonzero, the *reflect* flag is set if `z`
    should be replaced with `1-z` using the reflection formula.

    Note that this function does not guarantee the error bound rigorously;
    a rigorous error bound, which also accounts for the radius of `z`,
    is computed a posteriori when evaluating the Stirling series.
    However, in practice, this function does estimate the bound
    very accurately.

    To obtain a remainder smaller than `2^{-b}`, we must choose an `r` such
    that, in the real case, `x + r > \beta b`, where
    `\beta > \log(2) / (2 \pi) \approx 0.11`.
    In practice, a slightly larger factor `\beta \approx 0.2` more closely
    balances `n` and `r`. A much larger `\beta` (e.g. `\beta = 1`) could be
    used to reduce the number of Bernoulli numbers that have to be
    precomputed, at the expense of slower repeated evaluation.

.. function :: void gamma_stirling_coeff(fmprb_t b, ulong k, long prec)

    Sets `b = B_{2k} / (2k (2k-1))`, rounded to *prec* bits.
    Assumes that the Bernoulli number has been precomputed.

.. function :: void gamma_stirling_eval_series_fmprb(fmprb_t s, const fmprb_t z, long n, long prec)

.. function :: void gamma_stirling_eval_series_fmpcb(fmpcb_t s, const fmpcb_t z, long n, long prec)

    Evaluates the Stirling series

    .. math ::

        \log \Gamma(z) - R(n,z) = \left(z-\frac{1}{2}\right)\log z - z +
              \frac{\ln {2 \pi}}{2} + \sum_{k=1}^{n-1} t_k

    where

    .. math ::

        t_k = \frac{B_{2k}}{2k(2k-1)z^{2k-1}}.

    The bound

    .. math ::

        |R(n,z)| \le \frac{t_n}{\cos(0.5 \arg(z))^{2n}}

    is included in the radius of the output.


Rising factorials
--------------------------------------------------------------------------------

.. function :: void gamma_rising_fmprb_ui_bsplit_simple(fmprb_t y, const fmprb_t x, ulong n, long prec)

.. function :: void gamma_rising_fmprb_ui_bsplit_eight(fmprb_t y, const fmprb_t x, ulong n, long prec)

.. function :: void gamma_rising_fmprb_ui_bsplit_rectangular(fmprb_t y, const fmprb_t x, ulong n, ulong step, long prec)

.. function :: void gamma_rising_fmprb_ui_bsplit(fmprb_t y, const fmprb_t x, ulong n, long prec)

.. function :: void gamma_rising_fmpcb_ui_bsplit_simple(fmpcb_t y, const fmpcb_t x, ulong n, long prec)

.. function :: void gamma_rising_fmpcb_ui_bsplit_eight(fmpcb_t y, const fmpcb_t x, ulong n, long prec)

.. function :: void gamma_rising_fmpcb_ui_bsplit_rectangular(fmpcb_t y, const fmpcb_t x, ulong n, ulong step, long prec)

.. function :: void gamma_rising_fmpcb_ui_bsplit(fmpcb_t y, const fmpcb_t x, ulong n, long prec)

    Sets `y` to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`,
    computed using binary splitting.

    The different versions of this function process the basecase differently.
    The *simple* version simply multiplies together several factors
    one after another.

    The *eight* version processes eight factors at a time using the formula

    .. math ::

        x(x+1)\cdots(x+7) = (28 + 98x + 63x^2 + 14x^3 + x^4)^2 - 16 (7+2x)^2,

    replacing 7 full-precision multiplications with 3 squarings,
    1 multiplication, and several linear operations ([CP2005]_, page 316).
    Empirically, if `x` is a full-precision number, this is about twice as
    fast as the *simple* version at high precision. Numerical stability is
    slightly worse.

    The *rectangular* version processes *step* factors at a time by
    expanding the polynomial `f(t) = t (t+1) (t+2) \cdots (t+\mathrm{step}-1)`
    and evaluating each factor `f(x + \mathrm{step} \, k)`
    using rectangular splitting. At very high precision, if `x` is a
    full-precision number, this asymptotically reduces the number of
    full-precision multiplications required.

    The function *gamma_rising_fmprb_ui_bsplit* automatically chooses
    an algorithm depending on the inputs.

.. function :: void gamma_rising_fmprb_ui_multipoint(fmprb_t f, const fmprb_t c, ulong n, long prec)

    Sets `y` to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`,
    computed using fast multipoint evaluation. This only requires
    `O(n^{1/2+\varepsilon})` multiplications, but has high overhead
    and poor numerical stability (adding `O(n)` guard bits to the input
    might be necessary to achieve full accuracy). It can be expected to
    be faster than the binary splitting algorithm if the input is a
    full-precision number, the precision is at least 100000 bits,
    and *n* is of the same order of magnitude as (perhaps slightly
    smaller than) the number of bits.


Rational arguments
--------------------------------------------------------------------------------

.. function:: void gamma_series_fmpq_hypgeom(fmprb_struct * res, const fmpq_t a, long len, long prec)

    Given a rational number `0 < a \le 1`, uses binary splitting to compute
    *len* coefficients in the Taylor series of `\Gamma(a+x)`, i.e. computes
    `\Gamma(a), \Gamma'(a) ... \Gamma^{(\mathrm{len}-1)}(a) / (\mathrm{len}-1)!`.
    In particular, with *len* = 1, this function computes `\Gamma(a)`
    efficiently for small rational *a*.

    The *len* = 1 case of this algorithm dates back to Brent [Bre1978]_,
    and the extension to higher derivatives was done by Karatsuba [Kar1998]_.
    Karatsuba's original algorithm is suboptimal for large *len*;
    we use the faster algorithm given without error bounds
    by Borwein, Bradley and Crandall [BBC2000]_.

    The algorithm consists of evaluating the finite part of

    .. math ::

        \Gamma(s) = \int_0^{\infty} e^{-t} t^{s-1} dt = 
        N^s \left( \sum_{k=0}^R \frac{(-1)^k N^k}{(k + s) k!} + S \right) + I

    where

    .. math ::

        S = \sum_{k=R+1}^{\infty} \frac{(-1)^k N^k}{(k + s) k!}

    and

    .. math ::

        I = \int_N^{\infty} e^{-t} t^{s-1} dt.

    This formula is valid for complex `s` with `\Re{s} > 0`.
    It is therefore also valid if `s` is a power series argument `s = a + x`
    where `\Re{a} > 0`, so doing the arithmetic with truncated power
    series gives us the derivatives.

    We now discuss choosing the parameters `R` and `N`, and bounding
    the error terms `S` and `I`. We assume that `0 < \Re{a} \le 1`, `N \ge 1`
    and `R \ge 2 N`. The coefficients of `I` are given by

    .. math ::

        I = \int_N^{\infty} e^{-t} t^{a+x-1} dt =
        \sum_{j=0}^{\infty} \frac{x^j}{j!}
        \int_N^{\infty} e^{-t} t^{a-1} \log^j t \; dt.

    As shown by Karatsuba, the integrals are bounded in absolute value by
    `2 e^{-N} \log^j N`. Thus, for a precision of `p` bits, `N` should be
    about `p \log 2`.

    Expanding the terms in `S` as geometric series gives

    .. math ::

        S = \sum_{j=0}^{\infty} \, x^j \,
        \sum_{k=R+1}^{\infty} \frac{(-1)^{k+j} N^k}{(k+a)^{j+1} k!}.

    By the assumption that `R \ge 2 N`, the sums are bounded by

    .. math ::

        \frac{N^R}{R^{j+1} R!} \left(\frac{1}{2} + \frac{1}{4} + \ldots\right) =
        \frac{N^R}{R^{j+1} R!} \le \frac{1}{R^{j+1}} N^R \left(\frac{e}{R}\right)^R.

    Let `R = cN` where `c` is to be determined. Expanding the
    logarithm of `N^R \left(\frac{e}{R}\right)^R` around `N = \infty`
    gives the approximate magnitude `(c - c \log c) N`. Setting this equal
    to `-p \log 2`, we find that we should take
    `c = 1/W(1/e) \approx 3.59112147666862` where `W(x)` is the
    Lambert W-function. (Karatsuba gives the incorrect value `c = 3`).

    We also estimate the working precision needed in the binary splitting
    stage (the binary splitting could be done with exact arithmetic, but
    this is unnecessarily costly). Assume that the sum is around unity in
    magnitude. The binary logarithm of term `k` is roughly
    `b(k) = k \log_2 N + k \log_2 e - k \log_2 k`. Since
    `b'(k) = \log_2 N - \log_2 k`, the largest term magnitude occurs
    roughly at `k = N`, so we need to increase the working precision
    by about `b(N) = N / \log 2` bits.

