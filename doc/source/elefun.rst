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

.. function:: void elefun_exp_fixed_precomp(fmpz_t y, fmpz_t yerr, fmpz_t exponent, const fmpz_t x, const fmpz_t xerr, long prec)

    Given a fixed-point ball with midpoint *x* and radius *xerr*,
    computes a fixed-point ball *y* and *yerr* and a corresponding *exponent*
    giving the value of the exponential function.

    To avoid overflow, *xerr* must be small (for accuracy, it should be
    much smaller than 1).

    We first use division with remainder to remove `n \log 2`
    from `x`, leaving `0 \le x < 1` (and setting *exponent* to `n`).
    The value of `\log 2` is computed with 1 ulp error,
    which adds `|n|` ulp to the error of *x*.

    We then write `x = x_0 + x_1 + x_2` where `x_0` and `x_1` hold
    the top bits of `x`, looking up `E_0 = e^{x_0}` and
    `E_1 = e^{x_1}` from a table and computing
    `E_2 = e^{x_2}` using the Taylor series. Let the corresponding (absolute value)
    errors be `T_0, T_1, T_2` ulp (`\epsilon = 2^{-p}`).
    Next, we compute `E_0 (E_1 E_2)` using fixed-point multiplications.
    To bound the error, we write:

    .. math ::

        Y_0 = E_2 + T_2 \epsilon

        Y_1 = (E_1 + T_1 \epsilon) Y_0 + \epsilon

        Y_2 = (E_0 + T_0 \epsilon) Y_1 + \epsilon

    Expanding the expression for `Y_2 - E_0 E_1 E_2` and using
    `T_0, T_1 \le 1`, `E_0 \le 3`, `E_1, E_2 \le 1.1`,
    `\epsilon^2 \le \epsilon`, we find that the error is
    bounded by `9 + 7 T_2` ulp.

    Finally, we add the propagated error. We have
    `e^{a+b} - e^a \le b e^b e^a`. We bound `b e^b` by `b + b^2 + b^3`
    if `b \le 1`, and crudely by `4^b` otherwise.

.. function:: int elefun_exp_precomp(fmprb_t z, const fmprb_t x, long prec, int minus_one)

    Returns nonzero and sets *z* to a ball containing `e^x`
    (respectively `e^x-1` if *minus_one* is set), if *prec* and the input
    are small enough to efficiently (and accurately) use the fixed-point code
    with precomputation. If the precision or the arguments are too large,
    returns zero without altering *z*.

