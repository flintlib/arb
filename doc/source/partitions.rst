**partitions.h** -- computation of the partition function
===============================================================================

This module implements the asymptotically fast algorithm
for evaluating the integer partition function `p(n)`
described in [Joh2012]_.
The idea is to evaluate a truncation of the Hardy-Ramanujan-Rademacher series
using tight precision estimates, and symbolically factoring the
occurring exponential sums.

An implementation based on floating-point arithmetic can
also be found in FLINT. That version is significantly faster for
small `n` (e.g. `n < 10^6`), but relies on some numerical subroutines
that have not been proved correct.

The implementation provided here uses ball arithmetic throughout to guarantee
a correct error bound for the numerical approximation of `p(n)`.
For large `n`, it is nearly as fast as the floating-point version in FLINT.

.. function:: void partitions_rademacher_bound(fmpr_t b, ulong n, ulong N)

    Sets `b` to an upper bound for

    .. math ::

        M(n,N) = \frac{44 \pi^2}{225 \sqrt 3} N^{-1/2}
                  + \frac{\pi \sqrt{2}}{75} \left( \frac{N}{n-1} \right)^{1/2}
                \sinh\left(\frac{\pi}{N} \sqrt{\frac{2n}{3}}\right).

    This formula gives an upper bound for the truncation error in the
    Hardy-Ramanujan-Rademacher formula when the series is taken up
    to the term `t(n,N)` inclusive.

.. function:: void partitions_hrr_sum_fmprb(fmprb_t x, ulong n, long N0, long N)

    Evaluates the partial sum `\sum_{k=N_0}^N t(n,k)` of the
    Hardy-Ramanujan-Rademacher series.

.. function:: void partitions_fmpz_ui(fmpz_t p, ulong n)

    Computes the partition function `p(n)` using the Hardy-Ramanujan-Rademacher
    formula. This function computes a numerical ball containing `p(n)`
    and verifies that the ball contains a unique integer.

