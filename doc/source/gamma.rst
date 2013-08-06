.. _gamma:

**gamma.h** -- support for the gamma function
===============================================================================

This module implements various algorithms for evaluating the
gamma function and related functions. The functions provided here are mainly
intended for internal use, though they may be useful to call directly in some
applications where the default algorithm choices are suboptimal.
Most applications should use the user-friendly functions
in the :ref:`fmprb <fmprb>` and :ref:`fmpcb <fmpcb>` modules (or for
power series, the functions in the
:ref:`fmprb_poly <fmprb-poly>` and :ref:`fmpcb_poly <fmpcb-poly>`
modules).


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

.. function:: void gamma_stirling_choose_param_fmprb(int * reflect, long * r, long * n, const fmprb_t x, int allow_reflection, int digamma, long prec)

.. function:: void gamma_stirling_choose_param_fmpcb(int * reflect, long * r, long * n, const fmprb_t x, int allow_reflection, int digamma, long prec)

    Compute parameters `r`, `n` such that the remainder in the Stirling
    series with `z = x+yi` approximately is bounded by `2^{-\mathrm{prec}}`.
    If *digamma* is nonzero, the calculation is done for the digamma
    function.

    The parameter `n` is the truncation point in the asymptotic
    Stirling series. If `|z|` is too small for the Stirling series
    to give sufficient accuracy directly, we first translate to `z + r`
    using the formula `\Gamma(z) = \Gamma(z+r) / 
    (z (z+1) (z+2) \cdots (z+r-1))`.

    If *allow_reflection* is nonzero, the *reflect* flag is set if `z`
    should be replaced with `1-z` using the reflection formula.

    This function uses double precision arithmetic internally,
    and does not guarantee the error bound rigorously;
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

.. function :: void gamma_stirling_coeff(fmprb_t b, ulong k, int digamma, long prec)

    Sets `b = B_{2k} / (2k (2k-1))`, rounded to *prec* bits, or if *digamma*
    is nonzero, sets `b = B_{2k} / (2k)`.

.. function :: void gamma_stirling_eval_fmprb(fmprb_t s, const fmprb_t z, long n, int digamma, long prec)

.. function :: void gamma_stirling_eval_fmpcb(fmpcb_t s, const fmpcb_t z, long n, int digamma, long prec)

    Evaluates the Stirling series

    .. math ::

        \log \Gamma(z) = \left(z-\frac{1}{2}\right)\log z - z +
              \frac{\ln {2 \pi}}{2}
                + \sum_{k=1}^{n-1}  \frac{B_{2k}}{2k(2k-1)z^{2k-1}}
              + R(n,z).

    If *digamma* is nonzero, the derivative of this series (i.e. the
    expansion for the digamma function) is evaluated.
    The error bound for the tail `R(n,z)` (computed via
    :func:`gamma_stirling_bound_fmprb` or
    :func:`gamma_stirling_bound_fmpcb`) is included in the output.

.. function :: void gamma_stirling_eval_fmprb_series(fmprb_ptr res, const fmprb_t z, long n, long len, long prec)

.. function :: void gamma_stirling_eval_fmpcb_series(fmpcb_ptr res, const fmpcb_t z, long n, long len, long prec)

    Evaluates the Stirling series of a power series argument `z + t`,
    writing a power series truncated to length *len* to the output *res*.
    The error bound (computed via
    :func:`gamma_stirling_bound_fmprb` or
    :func:`gamma_stirling_bound_fmpcb`) is included in the output.

.. function :: void gamma_stirling_bound_phase(fmpr_t bound, const fmpcb_t z, long prec)

    Sets *bound* to an upper bound for the phase factor
    `b = 1/\cos(\operatorname{arg}(z)/2)` which appears in the error bound
    for the Stirling series. By trigonometric identities, assuming
    that `z = x+yi`, we have `b = \sqrt{1 + t^2}` where

    .. math ::

        t = \frac{y}{\sqrt{x^2 + y^2} + x} = \frac{\sqrt{x^2 + y^2} - x}{y}

    We bound `x` and `y` such that `|\operatorname{arg}(x+yi)|`
    is maximized, and then evaluate `t` with the choice of square root
    expression that avoids cancellation, using directional rounding throughout.

.. function :: void gamma_stirling_bound_fmprb(fmpr_struct * err, const fmprb_t z, long k0, long knum, long n)

.. function :: void gamma_stirling_bound_fmpcb(fmpr_struct * err, const fmpcb_t z, long k0, long knum, long n)

    Computes bounds for the truncation error in the Stirling series
    when summed up to term `n - 1` inclusive. An exact expression for the
    truncation error is given (see [Olv1997]_ pp. 293-295) by

    .. math ::

        R_n(z) = \int_0^{\infty} \frac{B_{2n} - {\tilde B}_{2n}(x)}{2n(x+z)^{2n}} dx.

    We optionally evaluate the bound for several terms in the
    Taylor series: considering `R_n(z+t) \in \mathbb{C}[[t]]`, we
    compute bounds for the coefficient of `t^k` for *knum* consecutive
    values of *k* starting with *k0*.
    Using the fact that the numerator of the integrand is bounded in
    absolute value by `2 |B_{2n}|`, and using the bound for `|x+z|`
    given by [Olv1997]_, we obtain

    .. math ::

        |[t^k] R_n(z+t)| \le 2 |B_{2n}|
            \frac{\Gamma(2n+k-1)}{\Gamma(k+1) \Gamma(2n+1)}
            \; |z| \; (b / |z|)^{2n+k}

    where `b` is the phase factor implemented by
    :func:`gamma_stirling_bound_phase`.


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
    full-precision multiplications required. If the *step* parameter
    is set to zero, a default value is used.

    The functions *gamma_rising_fmprb_ui_bsplit* and
    *gamma_rising_fmpcb_ui_bsplit* automatically choose
    an algorithm depending on the inputs.

.. function :: void gamma_rising_fmprb_ui_delta(fmprb_t y, const fmprb_t x, ulong n, ulong m, long prec)

.. function :: void gamma_rising_fmpcb_ui_delta(fmpcb_t y, const fmpcb_t x, ulong n, ulong m, long prec)

    Sets `y` to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`,
    computed as a product of partial products
    `(x+k)(x+k+1)\cdots(x+k+m-1)`. Each partial product is obtained
    from the previous by using a precomputed table of powers of `x` to
    evaluate the difference

    .. math ::

        \Delta_m(x,k) = (x+k+m)_{(m)} - (x+k)_{(m)}.

    The instance `m = 4` of this algorithm was used by Smith ([Smi2001]_),
    but we generalize it to a variable `m` which can be chosen nearly
    optimally depending on the precision and `n`.

    The polynomials `\Delta_m(x,k) \in \mathbb{Z}[k][x]` are generated dynamically.
    Expanding the rising factorials, applying the binomial theorem
    a couple of times, and doing several rearrangements of the sums, we
    find the closed form

    .. math ::

        \Delta_m(x,k) = \sum_{v=0}^{m-1} x^v \sum_{i=0}^{m-v-1} k^i C_m(v,i),

    where

    .. math ::

        C_m(v,i) = \sum_{j=i+1}^{m-v} m^{j-i} \left[{m \atop v+j}\right] {{v+j} \choose v} {j \choose i}

    in which the square bracket denotes an unsigned Stirling number
    of the first kind.


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

.. function :: void gamma_harmonic_sum_fmprb_ui_bsplit_simple(fmprb_t y, const fmprb_t x, ulong n, long prec)

.. function :: void gamma_harmonic_sum_fmprb_ui_bsplit_rectangular(fmprb_t y, const fmprb_t x, ulong n, ulong step, long prec)

.. function :: void gamma_harmonic_sum_fmprb_ui_bsplit(fmprb_t y, const fmprb_t x, ulong n, long prec)

.. function :: void gamma_harmonic_sum_fmpcb_ui_bsplit_simple(fmpcb_t y, const fmpcb_t x, ulong n, long prec)

.. function :: void gamma_harmonic_sum_fmpcb_ui_bsplit_rectangular(fmpcb_t y, const fmpcb_t x, ulong n, ulong step, long prec)

.. function :: void gamma_harmonic_sum_fmpcb_ui_bsplit(fmpcb_t y, const fmpcb_t x, ulong n, long prec)

    Sets `y` to the harmonic sum `1/x + 1/(x+1) + \ldots + 1/(x+n-1)`,
    computed using division-avoiding binary splitting.

    The *rectangular* version processes *step* terms at a time in
    analogy with the rising factorial algorithm.
    Letting `f(t) = t (t+1) (t+2) \cdots (t+n-1)`, we have
    `1/x + \ldots + 1/(x+n-1) = f'(x) / f(x)`.
    If the *step* parameter is set to zero, a default value is used.

    The functions *gamma_harmonic_sum_fmprb_ui_bsplit* and
    *gamma_harmonic_sum_fmpcb_ui_bsplit* automatically choose
    an algorithm depending on the inputs.


Rational arguments
--------------------------------------------------------------------------------

.. function:: void gamma_small_frac(fmprb_t y, unsigned int p, unsigned int q, long prec)

    Efficiently evaluates `y = \Gamma(p/q)` where `p/q` (assumed to be reduced)
    is one of `1, 1/2, 1/3, 2/3, 1/4, 3/4, 1/6, 5/6`.

    The cases `\Gamma(1) = 1` and `\Gamma(1/2) = \sqrt \pi` are trivial.
    We reduce all remaining cases to `\Gamma(1/3)` or `\Gamma(1/4)`
    using the following relations:

    .. math ::

        \Gamma(2/3) = \frac{2 \pi}{3^{1/2} \Gamma(1/3)}, \quad \quad
        \Gamma(3/4) = \frac{2^{1/2} \pi}{\Gamma(1/4)},

        \Gamma(1/6) = \frac{\Gamma(1/3)^2}{(\pi/3)^{1/2} 2^{1/3}}, \quad \quad
        \Gamma(5/6) = \frac{2 \pi (\pi/3)^{1/2} 2^{1/3}}{\Gamma(1/3)^2}.

    The values of `\Gamma(1/3)` and `\Gamma(1/4)` are cached for fast
    repeated evaluation. We compute them rapidly to high precision using

    .. math ::

        \Gamma(1/3) = \left( \frac{12 \pi^4}{\sqrt{10}}
            \sum_{k=0}^{\infty}
            \frac{(6k)!(-1)^k}{(k!)^3 (3k)! 3^k 160^{3k}} \right)^{1/6}, \quad \quad
        \Gamma(1/4) = \sqrt{\frac{(2\pi)^{3/2}}{\operatorname{agm}(1, \sqrt 2)}}.

    An alternative formula which could be used for `\Gamma(1/3)` is

    .. math ::

        \Gamma(1/3) = \frac{2^{4/9} \pi^{2/3}}{3^{1/12} \left( \operatorname{agm}\left(1,\frac{1}{2} \sqrt{2+\sqrt{3}}\right)\right)^{1/3}},

    but this appears to be slightly slower in practice.

.. function:: void gamma_series_fmpq_hypgeom(fmprb_ptr res, const fmpq_t a, long len, long prec)

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

