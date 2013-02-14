**gamma.h** -- support for the gamma function
===============================================================================

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


