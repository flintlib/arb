**elefun.h** -- support for evaluation of elementary functions
===============================================================================

The exponential function
--------------------------------------------------------------------------------

.. function:: long elefun_exp_taylor_bound(long mag, long prec)

    Returns *n* such that
    `\left|\sum_{k=n}^{\infty} x^k / k!\right| \le 2^{-\mathrm{prec}}`,
    assuming `|x| \le 2^{\mathrm{mag}} \le 1/4`.

.. function:: void elefun_exp_fixed_taylor_horner_precomp(fmpz_t y, fmpz_t h, const fmpz_t x, long n, long p)

    Letting `\epsilon = 2^{-p}`, computes `(y, h)` such that
    `|y \epsilon - \sum_{k=0}^{n-1} (x  \epsilon)^k / k!| \le h \epsilon`
    using Horner's rule and a precomputed table of
    inverse factorials. We assume that *n* and *p* are small enough
    and that `|x \epsilon| < 1`.

    The coefficients `c_k \approx 1/k!` are precomputed using truncating divisions,
    with error at most `2^{-q}` for some `q \ge p`. Rounding to `p` bits is
    then done using truncation. Since repeated truncation preserves
    correctness of rounding, the coefficients have error at most `\epsilon`.

    We now evaluate
    `s_0 = c_{n-1}, s_1 = R(x \epsilon s_0) + c_{n-2}, \ldots s_{n-1} = R(x \epsilon s_{n-2}) + c_0`.
    where `R(t)` denotes fixed-point truncation which adds an error
    of at most `\epsilon`. The error of `s_{i+1}` exceeds that of that
    of `s_i` by at most `2\epsilon`. Using induction, we obtain
    the bound `h = 2n-1`.

