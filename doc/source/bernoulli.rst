**bernoulli.h** -- support for Bernoulli numbers
===============================================================================

Generation of Bernoulli numbers
--------------------------------------------------------------------------------

.. type:: bernoulli_rev_t

.. function:: bernoulli_rev_init(bernrev_iter_t iter, ulong n)

.. function:: bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernrev_iter_t iter)

.. function:: bernoulli_rev_clear(bernrev_iter_t iter)

    Generates a range of even-indexed Bernoulli numbers in reverse order,
    i.e. generates `B_n, B_{n-2}, B_{n-4}, \ldots, B_0`.
    The exact, minimal numerators and denominators are produced.

    The Bernoulli numbers are computed by direct summation of the zeta series.
    This is made fast by storing a table of powers (as done by Bloemen et al.
    http://remcobloemen.nl/2009/11/even-faster-zeta-calculation.html).
    As an optimization, we only include the odd powers, and use
    fixed-point arithmetic.

    The reverse iteration order is preferred for performance reasons,
    as the powers can be updated using multiplications instead of divisions,
    and we avoid having to periodically recompute terms to higher precision.
    To generate Bernoulli numbers in the forward direction without having
    to store all of them, one can split the desired range into smaller
    blocks and compute each block with a single reverse pass.


Caching
-------------------------------------------------------------------------------

.. var:: long bernoulli_cache_num;

.. var:: fmpq * bernoulli_cache;

    Global cache of Bernoulli numbers.

.. function:: void bernoulli_cache_compute(long n)

    Makes sure that the Bernoulli numbers up to at least `B_{n-1}` are cached
    globally. Warning: this function is not currently threadsafe.


Bounding
-------------------------------------------------------------------------------

.. function:: long bernoulli_bound_2exp_si(ulong n)

    Returns an integer `b` such that `|B_n| \le 2^b`. Uses a lookup table
    for small `n`, and for larger `n` uses the inequality
    `|B_n| < 4 n! / (2 \pi)^n < 4 (n+1)^{n+1} e^{-n} / (2 \pi)^n`.
    Uses integer arithmetic throughout, with the bound for the logarithm
    being looked up from a table. If `|B_n| = 0`, returns *LONG_MIN*.
    Otherwise, the returned exponent `b` is never more than one percent
    larger than the true magnitude.

    This function is intended for use when `n` small enough that one might
    comfortably compute `B_n` exactly. It aborts if `n` is so large that
    internal overflow occurs.

