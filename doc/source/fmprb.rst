**fmprb.h** -- real numbers represented as floating-point balls
===============================================================================

An *fmprb_t* represents a ball over the real numbers,
that is, an interval `[m-r, m+r]` where the midpoint `m` and the
radius `r` are (extended) real numbers and `r` is nonnegative.
The result of an (approximate) operation done on *fmprb_t* variables
is a ball which contains the result of the (mathematically exact) operation
applied to any choice of points in the input balls.
In general, the output ball is not the smallest possible.

The precision parameter passed to each function roughly indicates the
precision to which calculations on the midpoint are carried out
(operations on the radius are always done using a fixed, small
precision.)

For arithmetic operations, the precision parameter currently simply
specifies the precision of the corresponding *fmpr* operation.
In the future, the arithmetic might be made faster by incorporating
sloppy rounding (typically equivalent to a loss of 1-2 bits of effective
working precision) when the result is known to be inexact (while still
propagating errors rigorously, of course).
Arithmetic operations done on exact input with exactly
representable output are always guaranteed to produce exact output.

For more complex operations, the precision parameter indicates a minimum
working precision (algorithms might allocate extra internal precision to
attempt to produce an output accurate to the requested number of bits,
especially when the required precision can be estimated easily, but this
is not generally required).

If the precision is increased and the inputs either are exact or are
computed with increased accuracy as well, the output should
converge proportionally, absent any bugs.
The general intended strategy for using ball arithmetic is to add a few
guard bits, and then repeat the calculation as necessary with an
exponentially increasing number of guard bits (Ziv's strategy) until the
result is exact
enough for one's purposes (typically the first attempt will be successful).
There are some caveats: in general, ball arithmetic only makes
sense for (Lipschitz) continuous functions, and 
trying to approximate functions close to singularities might result in
slow convergence, or failure to converge.

**Warning**: some methods for transcendental functions and constants
currently perform the error propagation in a non-rigorous way, due to the
implementation being incomplete (in some cases, a
rigorous error bound for the algorithm or function might not
be known at all).
This should be indicated in the documentation for each function.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmprb_struct

.. type:: fmprb_t

    An *fmprb_struct* consists of a pair of *fmpr_struct*:s.
    An *fmprb_t* is defined as an array of length one of type
    *fmprb_struct*, permitting an *fmprb_t* to be passed by
    reference.

.. macro:: FMPRB_RAD_PREC

    The precision used for operations on the radius. This is small
    enough to fit in a single word, currently 30 bits.

.. macro:: fmprb_midref(x)

    Macro returning a pointer to the midpoint of *x* as an *fmpr_t*.

.. macro:: fmprb_radref(x)

    Macro returning a pointer to the radius of *x* as an *fmpr_t*.


Memory management
-------------------------------------------------------------------------------

.. function:: void fmprb_init(fmprb_t x)

    Initializes the variable *x* for use. Its midpoint and radius are both
    set to zero.

.. function:: void fmprb_clear(fmprb_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: fmprb_struct * _fmprb_vec_init(long n)

    Returns a pointer to an array of *n* initialized *fmprb_struct*:s.

.. function:: void _fmprb_vec_clear(fmprb_struct * v, long n)

    Clears an array of *n* initialized *fmprb_struct*:s.


Basic manipulation
-------------------------------------------------------------------------------

.. function:: int fmprb_is_exact(const fmprb_t x)

    Returns nonzero iff the radius of *x* is zero.

.. function:: int fmprb_equal(const fmprb_t x, const fmprb_t y)

    Returns nonzero iff *x* and *y* are equal as balls, i.e. have both the
    same midpoint and radius.

.. function:: void fmprb_zero(fmprb_t x)

    Sets *x* to zero.

.. function:: int fmprb_is_zero(const fmprb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both zero.

.. function:: void fmprb_set(fmprb_t y, const fmprb_t x)

    Sets *y* to a copy of *x*.

.. function:: void fmprb_set_round(fmprb_t y, const fmprb_t x, long prec)

    Sets *y* to a copy of *x*, rounded to *prec* bits.

.. function:: void fmprb_neg(fmprb_t y, const fmprb_t x)

    Sets *y* to the negation of *x*.

.. function:: void fmprb_abs(fmprb_t y, const fmprb_t x)

    Sets *y* to the absolute value of *x*. No attempt is made to improve the
    interval represented by *x* if it contains zero.

.. function:: void fmprb_set_fmpr(fmprb_t y, const fmpr_t x)

.. function:: void fmprb_set_si(fmprb_t y, long x)

.. function:: void fmprb_set_ui(fmprb_t y, ulong x)

.. function:: void fmprb_set_fmpz(fmprb_t y, const fmpz_t x)

    Sets *y* exactly to *x*.

.. function:: void fmprb_set_fmpq(fmprb_t y, const fmpq_t x, long prec)

    Sets *y* to the rational number *x*, rounded to *prec* bits.

.. function:: void fmprb_one(fmprb_t x)

    Sets *x* to the exact integer 1.

.. function:: int fmprb_is_one(const fmprb_t x)

    Returns nonzero iff *x* is exactly 1.


Input and output
-------------------------------------------------------------------------------

.. function:: void fmprb_print(const fmprb_t x)

    Prints the internal representation of *x*.

.. function:: void fmprb_printd(const fmprb_t x, long digits)

    Prints *x* in decimal. The printed value of the radius is not adjusted
    to compensate for the fact that the binary-to-decimal conversion
    of both the midpoint and the radius introduces additional error.


Random number generation
-------------------------------------------------------------------------------

.. function:: void fmprb_randtest(fmprb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random ball. The midpoint and radius will both be finite.

.. function:: void fmprb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const fmprb_t x, long bits)

    Sets *q* to a random rational number from the interval represented by *x*.
    A denominator is chosen by multiplying the binary denominator of *x*
    by a random integer up to *bits* bits.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the exponent of the midpoint or radius is so large or small
    that representing the endpoints as exact rational numbers would
    cause overflows.


Precision and comparisons
-------------------------------------------------------------------------------

.. function:: void fmprb_add_error_fmpr(fmprb_t x, const fmpr_t err)

    Adds *err*, which is assumed to be nonnegative, to the radius of *x*.

.. function:: void fmprb_add_error_2exp_si(fmprb_t x, long e)

    Adds `2^e` to the radius of *x*.

.. function:: void fmprb_add_error(fmprb_t x, const fmprb_t err)

    Adds the supremum of *err*, which is assumed to be nonnegative, to the
    radius of *x*.

.. function:: int fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y)

.. function:: int fmprb_contains_fmpq(const fmprb_t x, const fmpq_t y)

.. function:: int fmprb_contains_fmpz(const fmprb_t x, const fmpz_t y)

.. function:: int fmprb_contains_mpfr(const fmprb_t x, const mpfr_t y)

.. function:: int fmprb_contains_zero(const fmprb_t x)

.. function:: int fmprb_contains(const fmprb_t x, const fmprb_t y)

    Returns nonzero iff the given number *y* is contained in the interval
    represented by *x*.

.. function:: int fmprb_overlaps(const fmprb_t x, const fmprb_t y)

    Returns nonzero iff *x* and *y* have some point in common.

.. function:: int fmprb_is_nonzero(const fmprb_t x)

    Returns nonzero iff zero is not contained in the interval represented
    by *x*.

.. function:: int fmprb_contains_negative(const fmprb_t x)

.. function:: int fmprb_contains_nonpositive(const fmprb_t x)

.. function:: int fmprb_contains_positive(const fmprb_t x)

.. function:: int fmprb_contains_nonnegative(const fmprb_t x)

    Returns nonzero iff there is any point *p* in the interval represented
    by *x* that is, respectively, `p < 0`, `p \le 0`, `p > 0`, `p \ge 0`.

.. function:: int fmprb_is_positive(const fmprb_t x)

.. function:: int fmprb_is_nonnegative(const fmprb_t x)

.. function:: int fmprb_is_negative(const fmprb_t x)

.. function:: int fmprb_is_nonpositive(const fmprb_t x)

    Returns nonzero iff all points *p* in the interval represented by *x*
    satisfy, respectively, `p > 0`, `p \ge 0`, `p < 0`, `p \le 0`.

.. void fmprb_get_abs_ubound_fmpr(fmpr_t u, const fmprb_t x, long prec)

    Sets *u* to the upper bound of the absolute value of *x*,
    rounded up to *prec* bits.

.. function:: void fmprb_get_abs_lbound_fmpr(fmpr_t u, const fmprb_t x, long prec)

    Sets *u* to the lower bound of the absolute value of *x*,
    rounded down to *prec* bits.

.. function:: void fmprb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const fmprb_t x)

    Computes the exact interval represented by *x*, in the form of an integer
    interval multiplied by a power of two, i.e. `x = [a, b] \times 2^{\mathrm{exp}}`.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the difference in magnitude between the midpoint and radius
    is so large that representing the endpoints exactly would cause overflows.

.. function:: int fmprb_get_unique_fmpz(fmpz_t z, const fmprb_t x)

    If *x* contains a unique integer, sets *z* to that value and returns
    nonzero. Otherwise (if *x* represents no integers or more than one integer),
    returns zero.

.. function:: long fmprb_rel_error_bits(const fmprb_t x)

    Returns the effective relative error of *x* measured in bits, defined as
    the difference between the position of the top bit in the radius
    and the top bit in the midpoint, plus one.
    The result is clamped between plus/minus *FMPR_PREC_EXACT*.

.. function:: long fmprb_rel_accuracy_bits(const fmprb_t x)

    Returns the effective relative accuracy of *x* measured in bits,
    equal to the negative of the return value from *fmprb_rel_error_bits*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmprb_add(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_add_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_add_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_add_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

.. function:: void fmprb_add_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)

    Sets `z = x + y`, rounded to *prec* bits. The precision can be
    *FMPR_PREC_EXACT* provided that the result fits in memory.

.. function:: void fmprb_sub(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_sub_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_sub_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_sub_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

    Sets `z = x - y`, rounded to *prec* bits. The precision can be
    *FMPR_PREC_EXACT* provided that the result fits in memory.

.. function:: void fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_mul_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

    Sets `z = x \times y`, rounded to *prec* bits. The precision can be
    *FMPR_PREC_EXACT* provided that the result fits in memory.

.. function:: void fmprb_mul_2exp_si(fmprb_t y, const fmprb_t x, long e)

    Sets *y* to *x* multiplied by `2^e`.

.. function:: void fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_div_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

.. function:: void fmprb_fmpz_div_fmpz(fmprb_t y, const fmpz_t num, const fmpz_t den, long prec)

.. function:: void fmprb_ui_div(fmprb_t z, ulong x, const fmprb_t y, long prec);

    Sets `z = x / y`, rounded to *prec* bits. If *y* contains zero, *z* is
    set to `0 \pm \infty`. Otherwise, error propagation uses the rule

    .. math ::
        \left| \frac{x}{y} - \frac{x+\xi_1 a}{y+\xi_2 b} \right| =
        \left|\frac{x \xi_2 b - y \xi_1 a}{y (y+\xi_2 b)}\right| \le
        \frac{|xb|+|ya|}{|y| (|y|-b)}

    where `-1 \le \xi_1, \xi_2 \le 1`, and
    where the triangle inequality has been applied to the numerator and
    the reverse triangle inequality has been applied to the denominator.

.. function:: void fmprb_div_2expm1_ui(fmprb_t y, const fmprb_t x, ulong n, long prec);

    Sets `y = x / (2^n - 1)`, rounded to *prec* bits.

.. function:: void fmprb_addmul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_addmul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_addmul_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_addmul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

    Sets `z = z + x \times y`, rounded to prec bits. The precision can be
    *FMPR_PREC_EXACT* provided that the result fits in memory.

.. function:: void fmprb_submul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_submul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_submul_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_submul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

    Sets `z = z - x \times y`, rounded to *prec* bits. The precision can be
    *FMPR_PREC_EXACT* provided that the result fits in memory.

.. function:: void fmprb_sqrt(fmprb_t z, const fmprb_t x, long prec)

.. function:: void fmprb_sqrt_ui(fmprb_t z, ulong x, long prec)

.. function:: void fmprb_sqrt_fmpz(fmprb_t z, const fmpz_t x, long prec)

    Sets *z* to the square root of *x*, rounded to *prec* bits.
    Error propagation is done using the following rule:
    assuming `m > r \ge 0`, the error is largest at `m - r`, and we have
    `\sqrt{m} - \sqrt{m-r} \le r / (2 \sqrt{m-r})`.

.. function:: void fmprb_sqrtpos(fmprb_t z, const fmprb_t x, long prec)

    Sets *z* to the square root of *x*, assuming that *x* represents a
    nonnegative number (i.e. discarding any negative numbers in the input
    interval), and producing an output interval not containing any
    negative numbers (unless the radius is infinite).

.. function:: void fmprb_pow_fmpz(fmprb_t y, const fmprb_t b, const fmpz_t e, long prec)

.. function:: void fmprb_pow_ui(fmprb_t y, const fmprb_t b, ulong e, long prec)

.. function:: void fmprb_ui_pow_ui(fmprb_t y, ulong b, ulong e, long prec)

.. function:: void fmprb_si_pow_ui(fmprb_t y, long b, ulong e, long prec)

    Sets `y = b^e` using binary exponentiation. Provided that *b* and *e*
    are small enough and the exponent is positive, the exact power can be
    computed using *FMPR_PREC_EXACT*.


Special functions
-------------------------------------------------------------------------------

Elementary functions
...............................................................................

.. function:: void fmprb_log(fmprb_t z, const fmprb_t x, long prec)

.. function:: void fmprb_log_ui(fmprb_t z, ulong x, long prec)

.. function:: void fmprb_log_fmpz(fmprb_t z, const fmpz_t x, long prec)

    Sets `z = \log(x)`. Error propagation is done using the following rule:
    assuming `m > r \ge 0`, the error is largest at `m - r`, and we have
    `\log(m) - \log(m-r) = \log(1 + r/(m-r))`. The last expression is
    calculated accurately for small radii via *fmpr_log1p*.
    An input containing zero currently raises an exception.

.. function:: void fmprb_exp(fmprb_t z, const fmprb_t x, long prec)

    Sets `z = \exp(x)`. Error propagation is done using the following rule:
    the error is largest at `m + r`, and we have
    `\exp(m+r) - \exp(m) = \exp(m) (\exp(r)-1)`.
    The last expression is calculated accurately for small radii
    via *fmpr_expm1*.

.. function:: void fmprb_sin(fmprb_t s, const fmprb_t x, long prec)

.. function:: void fmprb_cos(fmprb_t c, const fmprb_t x, long prec)

.. function:: void fmprb_sin_cos(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)

    Sets `s = \sin x`, `c = \cos x`. Error propagation uses the rule
    `|\sin(m \pm r) - \sin(m)| \le r` (this could be tightened to
    `\min(r,2)`).

.. function:: void fmprb_atan(fmprb_t z, const fmprb_t x, long prec)

    Sets `z = \tan^{-1} x`. Letting `d = \max(0, |m| - r)`,
    the propagated error is bounded by `r / (1 + d^2)`
    (this could be tightened).

.. function:: void fmprb_sinh(fmprb_t s, const fmprb_t x, long prec)

.. function:: void fmprb_cosh(fmprb_t c, const fmprb_t x, long prec)

.. function:: void fmprb_sinh_cosh(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)

    Sets `s = \sinh x`, `c = \cosh x`. Error propagation uses
    the derivatives of the functions.

Factorials and other integer functions
...............................................................................

.. function:: void fmprb_fac_ui(fmprb_t x, ulong n, long prec)

    Sets *x* to *n* factorial, computed using binary splitting. Provided that
    *n* is small enough, the exact factorial can be computed using
    *FMPR_PREC_EXACT*.

.. function:: void fmprb_rfac_ui_bsplit(fmprb_t y, const fmprb_t x, ulong n, long prec)

    Sets *x* to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`,
    computed using binary splitting.
    The basecase processes eight factors at a time using the formula

    .. math ::

        x(x+1)\cdots(x+7) = (28 + 98x + 63x^2 + 14x^3 + x^4)^2 - 16 (7+2x)^2,

    replacing 7 full-precision multiplications with 4 squarings and
    1 multiplication ([CP2005]_, page 316).
    Empirically, this is about twice as fast at high precision.
    Numerical stability is slightly worse.

.. function:: void fmprb_rfac_ui_multipoint(fmprb_t y, const fmprb_t x, ulong n, long prec)

    Sets *x* to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`,
    computed using fast multipoint evaluation. This only requires
    `O(n^{1/2+\varepsilon})` multiplications, but has high overhead
    and poor numerical stability (adding `O(n)` guard bits to the input
    might be necessary to achieve full accuracy). It can be expected to
    be faster than the binary splitting algorithm if the input is a
    full-precision number, the precision is at least 100000 bits,
    and *n* is of the same order of magnitude as (perhaps slightly
    smaller than) the number of bits.

.. function:: void fmprb_bin_ui(fmprb_t x, const fmprb_t n, ulong k, long prec)

.. function:: void fmprb_bin_uiui(fmprb_t x, ulong n, ulong k, long prec)

    Sets *x* to the binomial coefficient `{n \choose k}`, computed using
    binary splitting. Provided that *n* and *k* are small enough, an exact
    binomial coefficient can be computed using *FMPR_PREC_EXACT*.

.. function:: void fmprb_fib_fmpz(fmprb_t f, const fmpz_t n, long prec)

.. function:: void fmprb_fib_ui(fmprb_t f, ulong n, long prec)

    Sets x to the Fibonacci number `F_n`. Uses the binary squaring
    algorithm described in [Tak2000]_.
    Provided that *n* is small enough, an exact Fibonacci number can be
    computed using *FMPR_PREC_EXACT*.

Constants
...............................................................................

.. function:: void fmprb_const_pi_chudnovsky(fmprb_t x, long prec)

    Sets *x* to `\pi`, computed using the Chudnovsky algorithm.
    Letting `A = 13591409`, `B =  545140134`, `C = 640320`,
    we have `\pi \approx 1 / s_N` where

    .. math ::

        s_N = 12 \sum_{k=0}^N \frac{(-1)^k (6k)! (A+Bk)}
            {(3k)! (k!)^3 C^{3k+3/2}}

    The implementation computes an approximation for the
    algebraic number `1/s_N` using binary splitting, bounding
    the rounding error automatically.
    The hypergeometric term ratio is asymptotically
    `R = C^3 / (2^6 \times 3^3) \approx 1.5 \times 10^{14}`, and in fact we have
    `|\pi - 1/s_N| < 1/R^N` (with a more detailed calculation, the truncation
    error could be bounded closer to `1/R^{N+1}`).

.. function:: void fmprb_const_pi(fmprb_t x, long prec)

    Sets *x* to `\pi`. The value is cached for repeated use.

.. function:: void fmprb_const_log_sqrt2pi(fmprb_t x, long prec)

    Sets *x* to `\log \sqrt{2 \pi}`. The value is cached for repeated use.

.. function:: void fmprb_const_euler_brent_mcmillan(fmprb_t res, long prec)

    Sets *x* to Euler's constant `\gamma`, computed using the second
    Bessel function formula of Brent and McMillan ([BM1980]_,  [MPFR2012]_).
    Brent and McMillan conjectured that the error depending
    on the internal parameter *n* is of order `O(e^{-8n})`. Brent has
    recently proved that this bound is correct, but without determining
    an explicit big-O factor [Bre2010]_.

.. function:: void fmprb_const_khinchin(fmprb_t res, long prec)

    Sets *res* to Khinchin's constant `K_0`, computed as

    .. math ::

        \log K_0 = \frac{1}{\log 2} \left[
        \sum_{k=2}^{N-1} \log \left(\frac{k-1}{k} \right) \log \left(\frac{k+1}{k} \right)
        + \sum_{n=1}^\infty 
        \frac {\zeta (2n,N)}{n} \sum_{k=1}^{2n-1} \frac{(-1)^{k+1}}{k}
        \right]

    where `N \ge 2` is a free parameter that can be used for tuning [BBC1997]_.
    If the infinite series is truncated after `n = M`, the remainder
    is smaller in absolute value than

    .. math ::

        \sum_{n=M+1}^{\infty} \zeta(2n, N) = 
        \sum_{n=M+1}^{\infty} \sum_{k=0}^{\infty} (k+N)^{-2n} \le
        \sum_{n=M+1}^{\infty} \left( N^{-2n} + \int_0^{\infty} (t+N)^{-2n} dt \right)

        = \sum_{n=M+1}^{\infty} \frac{1}{N^{2n}} \left(1 + \frac{N}{2n-1}\right)
        \le \sum_{n=M+1}^{\infty} \frac{N+1}{N^{2n}} = \frac{1}{N^{2M} (N-1)}
        \le \frac{1}{N^{2M}}.

    Thus, for an error of at most `2^{-p}` in the series,
    it is sufficient to choose `M \ge p / (2 \log_2 N)`.

Riemann zeta function
...............................................................................

.. function:: void fmprb_const_zeta3_bsplit(fmprb_t x, long prec)

    Sets *x* to Apery's constant `\zeta(3)`, computed by applying binary
    splitting to a hypergeometric series.

.. function:: void fmprb_zeta_ui_asymp(fmprb_t z, ulong s, long prec)

    Assuming `s \ge 2`, approximates `\zeta(s)` by `1 + 2^{-s}` along with
    a correct error bound. We use the following bounds: for `s > b`,
    `\zeta(s) - 1 < 2^{-b}`, and generally,
    `\zeta(s) - (1 + 2^{-s}) < 2^{2-\lfloor 3 s/2 \rfloor}`.

.. function:: void fmprb_zeta_ui_euler_product(fmprb_t z, ulong s, long prec)

    Computes `\zeta(s)` using the Euler product. This is fast only if *s*
    is large compared to the precision.

    Writing `P(a,b) = \prod_{a \le p \le b} (1 - p^{-s})`, we have
    `1/\zeta(s) = P(a,M) P(M+1,\infty)`.

    To bound the error caused by truncating
    the product at `M`, we write `P(M+1,\infty) = 1 - \epsilon(s,M)`.
    Since `0 < P(a,M) \le 1`, the absolute error for `\zeta(s)` is
    bounded by `\epsilon(s,M)`.

    According to the analysis in [Fil1992]_, it holds for all `s \ge 6` and `M \ge 1`
    that `1/P(M+1,\infty) - 1 \le f(s,M) \equiv 2 M^{1-s} / (s/2 - 1)`.
    Thus, we have `1/(1-\epsilon(s,M)) - 1 \le f(s,M)`, and expanding
    the geometric series allows us to conclude that
    `\epsilon(M) \le f(s,M)`.

.. function:: void fmprb_zeta_ui_bernoulli(fmprb_t x, ulong n, long prec)

    Computes `\zeta(n)` for even *n* via the corresponding Bernoulli number,
    which is generated using FLINT.

.. function:: void fmprb_zeta_ui_vec_borwein(fmprb_struct * z, ulong start, long num, ulong step, long prec)

    Evaluates `\zeta(s)` at `\mathrm{num}` consecutive integers *s* beginning
    with *start* and proceeding in increments of *step*.
    Uses Borwein's formula ([Bor2000]_, [GS2003]_),
    implemented to support fast multi-evaluation
    (but also works well for a single *s*).

    Requires `\mathrm{start} \ge 2`. For efficiency, the largest *s*
    should be at most about as
    large as *prec*. Arguments approaching *LONG_MAX* will cause
    overflows.
    One should therefore only use this function for *s* up to about *prec*, and
    then switch to the Euler product.

    The algorithm for single *s* is basically identical to the one used in MPFR
    (see [MPFR2012]_ for a detailed description).
    In particular, we evaluate the sum backwards to avoid storing more than one
    `d_k` coefficient, and use integer arithmetic throughout since it
    is convenient and the terms turn out to be slightly larger than
    `2^\mathrm{prec}`.
    The only numerical error in the main loop comes from the division by `k^s`,
    which adds less than 1 unit of error per term.
    For fast multi-evaluation, we repeatedly divide by `k^{\mathrm{step}}`.
    Each division reduces the input error and adds at most 1 unit of
    additional rounding error, so by induction, the error per term
    is always smaller than 2 units.

.. function:: void fmprb_zeta_ui_bsplit(fmprb_t x, ulong s, long prec)

    Computes `\zeta(s)` for arbitrary `s \ge 2` using a binary splitting
    implementation of Borwein's formula. The algorithm has quasilinear
    complexity with respect to the precision.

.. function:: void fmprb_zeta_ui(fmprb_t x, ulong s, long prec)

    Computes `\zeta(s)` for nonnegative integer `s \ne 1`, automatically
    choosing an appropriate algorithm.

.. function:: void fmprb_zeta_ui_vec(fmprb_struct * x, ulong start, long num, long prec)

.. function:: void fmprb_zeta_ui_vec_even(fmprb_struct * x, ulong start, long num, long prec)

.. function:: void fmprb_zeta_ui_vec_odd(fmprb_struct * x, ulong start, long num, long prec)

    Computes `\zeta(s)` at num consecutive integers (respectively num
    even or num odd integers) beginning with `s = \mathrm{start} \ge 2`,
    automatically choosing an appropriate algorithm.

Gamma function
...............................................................................

.. function:: void fmprb_gamma_fmpq_karatsuba(fmprb_struct * v, const fmpq_t a, long num, long prec)

    Uses Karatsuba's algorithm [Kar1998]_ to compute num coefficients in the
    Taylor series of `\Gamma(a+x)` for rational `0 < a \le 1`, i.e.
    computes `\Gamma(a), \Gamma'(a) ... \Gamma^{(\mathrm{num}-1)}(a) / (\mathrm{num}-1)!`
    This algorithm is most efficient at high precision, for num much smaller
    than the number of bits, and with small denominators of *a*.
    In particular, with num = 1, this algorithm computes `\Gamma(a)`
    efficiently for small rational *a*.

    Let `s = \max(2, \mathrm{num}-1)`. With parameters `r` and `n` chosen
    such that `r \ge n` and `n \ge 2 s \log 2 s`, Karatsuba shows that

    .. math ::

        \Gamma^{(j)}(a) = \sum_{k=0}^r \frac{(-1)^k}{k!}
        \frac{n^{k+a}}{k+a}
        \sum_{m=0}^j (-1)^m \frac{j! \, \log^{j-m} n}{(j-m)! (k+a)^m}
        + \theta_j

    where

    .. math ::

        |\theta_j| \le \frac{5}{3} e^{-n} \log^s n +
        \left(\frac{e}{r+2}\right)^{r+2} (1 + n^{r+2} \log^s n).

    We choose the parameters `n` and `r` heuristically to be nearly optimal,
    and then evaluate the above formula to bound `\theta_j` rigorously.

    Karatsuba claims that choosing `r \ge 3n` gives
    `|\theta_j| \le 2^{-n-1}`. This is, unfortunately, incorrect.
    Setting `r = n \alpha` and expanding the error term around `n = \infty`,
    one finds that `\alpha` asymptotically should be
    `1/W(1/e) \approx 3.59112147666862` where `W(x)` is the Lambert W-function.
    We also optimize the selection of `n` by choosing
    `n \approx b \log 2` where `b` is the desired number of bits, rather
    than `n \approx b`, and round `n` so that it has a short binary expansion
    (this gives smaller numbers in the binary splitting stage).

    Finally, if `s` is small, we perform binary splitting to a working
    precision of about `2.2` times the target precision rather than exactly.
    This factor was tested to give full accuracy up to at least one million
    digits when `s \approx 1`. A more careful analysis should
    be done here so that a working precision is selected which always is
    sufficient and also nearly optimal.

.. function:: void fmprb_gamma_log(fmprb_t y, const fmprb_t x, long prec)

    Sets `y = \log \Gamma(x)`, assuming that `x > 0`.

    For large `x`, uses Stirling's expansion

    .. math ::

        \log \Gamma(x) = \left(x-\frac{1}{2}\right)\log x - x +
        \frac{\ln {2 \pi}}{2} + \sum_{k=1}^{n-1} t_k + R(n,x)

    where

    .. math ::

        t_k = \frac{B_{2k}}{2k(2k-1)x^{2k-1}}

    and `|R(n,x)| < t_n`.

    If `x` is too small for the asymptotic expansion to give sufficient
    accuracy directly, we translate to `x + r`
    using the formula `\log \Gamma(x) = \log \Gamma(x+r) -
    \log(x (x+1) (x+2) \cdots (x+r-1))`.

    To obtain a remainder smaller than `2^{-b}`, we must choose an `r` such
    that `x + r > \beta b`, where `\beta > \log(2) / (2 \pi) \approx 0.11`.
    We use a slightly larger factor `\beta \approx 0.2` to more closely
    balance `n` and `r`. A much larger `\beta` (e.g. `\beta = 1`) could be
    used to reduce the number of Bernoulli numbers that have to be
    precomputed, at the expense of slower repeated evaluation.

