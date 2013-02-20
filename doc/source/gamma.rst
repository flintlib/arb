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


