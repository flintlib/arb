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

.. function:: void fmprb_set_fmpz_2exp(fmprb_t x, const fmpz_t y, const fmpz_t exp)

    Sets *x* to *y* multiplied by 2 raised to the power *exp*.

.. function:: void fmprb_set_round_fmpz_2exp(fmprb_t y, const fmpz_t x, const fmpz_t exp, long prec)

    Sets *x* to *y* multiplied by 2 raised to the power *exp*, rounding
    the result to *prec* bits.


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

.. function:: void fmprb_randtest_exact(fmprb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random number with zero radius.

.. function:: void fmprb_randtest_precise(fmprb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random number with radius at most `2^{-\mathrm{prec}}`
    the magnitude of the midpoint.

.. function:: void fmprb_randtest_wide(fmprb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random number with midpoint and radius chosen independently,
    possibly giving a very large interval.

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

.. function:: void fmprb_set_interval_fmpr(fmprb_t x, const fmpr_t a, const fmpr_t b, long prec)

    Sets *x* to a ball containing the interval `[a, b]`. We
    require that `a \le b`.

.. function:: void fmprb_union(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

    Sets *z* to a ball containing both *x* and *y*.


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

.. function:: void fmprb_ui_div(fmprb_t z, ulong x, const fmprb_t y, long prec)

    Sets `z = x / y`, rounded to *prec* bits. If *y* contains zero, *z* is
    set to `0 \pm \infty`. Otherwise, error propagation uses the rule

    .. math ::
        \left| \frac{x}{y} - \frac{x+\xi_1 a}{y+\xi_2 b} \right| =
        \left|\frac{x \xi_2 b - y \xi_1 a}{y (y+\xi_2 b)}\right| \le
        \frac{|xb|+|ya|}{|y| (|y|-b)}

    where `-1 \le \xi_1, \xi_2 \le 1`, and
    where the triangle inequality has been applied to the numerator and
    the reverse triangle inequality has been applied to the denominator.

.. function:: void fmprb_div_2expm1_ui(fmprb_t y, const fmprb_t x, ulong n, long prec)

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

Powers and roots
-------------------------------------------------------------------------------

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

.. function:: void fmprb_root(fmprb_t z, const fmprb_t x, ulong k, long prec)

    Sets *z* to the *k*-th root of *x*, rounded to *prec* bits.
    Warning: this function is only fast for small fixed *k*. For large *k*,
    it is better to use the exponential function.

.. function:: void fmprb_pow_fmpz(fmprb_t y, const fmprb_t b, const fmpz_t e, long prec)

.. function:: void fmprb_pow_ui(fmprb_t y, const fmprb_t b, ulong e, long prec)

.. function:: void fmprb_ui_pow_ui(fmprb_t y, ulong b, ulong e, long prec)

.. function:: void fmprb_si_pow_ui(fmprb_t y, long b, ulong e, long prec)

    Sets `y = b^e` using binary exponentiation. Provided that *b* and *e*
    are small enough and the exponent is positive, the exact power can be
    computed using *FMPR_PREC_EXACT*.

.. function:: void fmprb_agm(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

    Sets *z* to the arithmetic-geometric mean of *x* and *y*.

Exponentials and logarithms
-------------------------------------------------------------------------------

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


Trigonometric functions
-------------------------------------------------------------------------------

.. function:: void fmprb_sin(fmprb_t s, const fmprb_t x, long prec)

.. function:: void fmprb_cos(fmprb_t c, const fmprb_t x, long prec)

.. function:: void fmprb_sin_cos(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)

    Sets `s = \sin x`, `c = \cos x`. Error propagation uses the rule
    `|\sin(m \pm r) - \sin(m)| \le r` (this could be tightened to
    `\min(r,2)`).

.. function:: void fmprb_sin_pi(fmprb_t s, const fmprb_t x, long prec)

.. function:: void fmprb_cos_pi(fmprb_t c, const fmprb_t x, long prec)

.. function:: void fmprb_sin_cos_pi(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)

    Sets `s = \sin \pi x`, `c = \cos \pi x`.

.. function:: void fmprb_sin_pi_fmpq(fmprb_t s, const fmpq_t x, long prec)

.. function:: void fmprb_cos_pi_fmpq(fmprb_t c, const fmpq_t x, long prec)

.. function:: void fmprb_sin_cos_pi_fmpq(fmprb_t s, fmprb_t c, const fmpq_t x, long prec)

    Sets `s = \sin \pi x`, `c = \cos \pi x` where `x` is a rational
    number (whose numerator and denominator are assumed to be reduced).
    We first use trigonometric symmetries to reduce the argument to the
    octant `[0, 1/4]`. Then we either multiply by a numerical approximation
    of `\pi` and evaluate the trigonometric function the usual way,
    or we use algebraic methods (*_fmprb_sin_pi_fmpq_algebraic* et al),
    depending on which is estimated to be faster.
    Since the argument has been reduced to the first octant, the
    first of these two methods gives full accuracy even if the original
    argument is close to some root other the origin.

.. function:: void _fmprb_sin_pi_fmpq_algebraic(fmprb_t s, ulong p, ulong q, long prec)

.. function:: void _fmprb_cos_pi_fmpq_algebraic(fmprb_t c, ulong p, ulong q, long prec)

.. function:: void _fmprb_sin_cos_pi_fmpq_algebraic(fmprb_t s, fmprb_t c, ulong p, ulong q, long prec)

    Uses algebraic methods to evaluate `s = \sin(p \pi / q)`,
    `c = \cos(p \pi / q)` where `0 \le 2p \le q` and `\gcd(p,q) = 1`.
    This is efficient if `q` has the form `2^r`, `3 \times 2^r` or `5 \times 2^r`,
    with `r \ge 0`, or if `q` is a moderately large integer and the precision
    is in the thousands of bits (otherwise simply evaluating
    the trigonometric function as a transcendental is cheaper).

    We use direct formulas if `1 \le q \le 6`.
    Otherwise, consider the cosine case (we shift the sine into a cosine,
    and for evaluating both functions simultaneously, we use the Pythagorean
    theorem `\sin x = \pm \sqrt{1-\cos^2 x}`, costing one extra square root).

    We first remove the largest power of two `2^r` dividing `q` by repeatedly
    doubling the angle (requiring the computation of `r` nested square roots).
    If `q = 2^r` or `q = 3 \times 2^r` or `q = 5 \times 2^r` this allows us to
    recurse all the way to the direct formulas, and we are done.

    Otherwise, having transformed `p, q` so that `q` is odd,
    we generate the minimal polynomial in `\mathbb{Z}[x]` of the
    algebraic number `\cos(p \pi / q)` and refine a low-precision
    value of the root to high accuracy using Newton iteration.

    This function assumes that `q` is small for correct operation.
    In particular, it assumes that `4p` does not overflow a limb.
    For efficiency, we also assume that `q / 2^r` is reasonably small
    (otherwise the minimal polynomial becomes impractically large, possibly
    exhausting the available memory).

.. function:: void fmprb_atan(fmprb_t z, const fmprb_t x, long prec)

    Sets `z = \tan^{-1} x`. Letting `d = \max(0, |m| - r)`,
    the propagated error is bounded by `r / (1 + d^2)`
    (this could be tightened).

.. function:: void fmprb_atan2(fmprb_t r, const fmprb_t b, const fmprb_t a, long prec)

    Sets *r* to an the argument (phase) of the complex number
    `a + bi`, with the branch cut discontinuity on `(-\infty,0]`.
    We define `\operatorname{atan2}(0,0) = 0`, and for `a < 0`,
    `\operatorname{atan2}(0,a) = \pi`.

Hyperbolic functions
-------------------------------------------------------------------------------

.. function:: void fmprb_sinh(fmprb_t s, const fmprb_t x, long prec)

.. function:: void fmprb_cosh(fmprb_t c, const fmprb_t x, long prec)

.. function:: void fmprb_sinh_cosh(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)

    Sets `s = \sinh x`, `c = \cosh x`. Error propagation uses
    the derivatives of the functions.


Factorials and other integer functions
-------------------------------------------------------------------------------

.. function:: void fmprb_fac_ui(fmprb_t x, ulong n, long prec)

    Sets *x* to *n* factorial, computed using binary splitting. Provided that
    *n* is small enough, the exact factorial can be computed using
    *FMPR_PREC_EXACT*.

.. function:: void fmprb_rising_fmprb_ui(fmprb_t y, const fmprb_t x, ulong n, long prec)

    Sets *x* to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`.

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
-------------------------------------------------------------------------------

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
-------------------------------------------------------------------------------

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
    implementation of Borwein's algorithm. This has quasilinear complexity
    with respect to the precision (assuming that `s` is fixed).
    We have

    .. math ::

        \zeta(s) = \frac{1}{d_n (1-2^{1-s})}
        \sum_{k=0}^{n-1} \frac{(-1)^k(d_n-d_k)}{(k+1)^s} + \gamma_n(s)

    where

    .. math ::

        d_k = n \sum_{i=0}^k \frac{(n+i-1)! 4^i}{(n-i)! (2i)!}.

    On the domain of interest, `|\gamma_n(s)| \le 3 / (3 + \sqrt 8)^n`.

    We write the summation as a system of first-order recurrences for
    `(s_k, d_k, t_k)` where `t_k = d_k - d_{k-1}`. This system is
    described by the matrix equation

    .. math ::

        \begin{pmatrix} s_{k+1} \\ d_{k+2} \\ t_{k+3} \end{pmatrix}
        =
        \begin{pmatrix}
        1 & (-1)^k (k+1)^{-s} & 0 \\
        0 & 1 & 1 \\
        0 & 0 & u(k)
        \end{pmatrix}
        \begin{pmatrix} s_k \\ d_{k+1} \\ t_{k+2} \end{pmatrix}.

    We derive the binary splitting scheme by considering a product
    of an arbitrary pair in the chain `M_{n-1} M_{n-2} \cdots M_1 M_0`.
    This gives

    .. math ::

        \begin{pmatrix}
        1 & A_L & B_L \\
        0 & 1 & C_L \\
        0 & 0 & D_L
        \end{pmatrix}
        \begin{pmatrix}
        1 & A_R & B_R \\
        0 & 1 & C_R \\
        0 & 0 & D_R
        \end{pmatrix} =
        \begin{pmatrix}
        1 & A_L+A_R & B_R+A_L C_R+B_L D_R \\
        0 & 1 & C_R+C_L D_R \\
        0 & 0 & D_L D_R
        \end{pmatrix}.

    The next step is to clear denominators. Instead of putting the
    whole matrix on a common denominator, we optimize by putting `C, D` on a
    denominator `Q_1` (the product of denominators of `u`) and `A, B` on
    a common denominator `Q_3 = Q_1 Q_2` (where `Q_2` is the product of
    `(k+1)^s` factors). This gives a small efficiency improvement. Thus,
    we have

    .. math ::

        \begin{pmatrix}
        1 & \dfrac{A_L}{Q_{3L}} & \dfrac{B_L}{Q_{3L}} \\[3ex]
        0 & 1 & \dfrac{C_L}{Q_{1L}} \\[3ex]
        0 & 0 & \dfrac{D_L}{Q_{1L}}
        \end{pmatrix}
        \begin{pmatrix}
        1 & \dfrac{A_R}{Q_{3R}} & \dfrac{B_R}{Q_{3R}} \\[3ex]
        0 & 1 & \dfrac{C_R}{Q_{1R}} \\[3ex]
        0 & 0 & \dfrac{D_R}{Q_{1R}}
        \end{pmatrix} =
        \begin{pmatrix}
        1 & \dfrac{Q_{3L} A_R + A_L Q_{3R}}{Q_{3L} Q_{3R}} & \dfrac{Q_{3L} B_R + A_L C_R Q_{2R} + B_L D_R Q_{2R}}{Q_{3L} Q_{3R}} \\[3ex]
        0 & 1 & \dfrac{Q_{1L} C_R + C_L D_R}{Q_{1L} Q_{1R}} \\[3ex]
        0 & 0 & \dfrac{D_L D_R}{Q_{1L} Q_{1R}}
        \end{pmatrix}.

    In the final matrix, we note that 
    `A / Q_3 = \sum_k (-1)^k (k+1)^{-s}`, and `C / Q_1 = d_n`.
    Thus `(1 / d_n) \sum_k (-1)^k (k+1)^{-s} (d_n - d_k)` is given by
    `A/Q_3 - (B/Q_3) / (C/Q_1) = (A C - B Q_1) / (Q_3 C)`.

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
-------------------------------------------------------------------------------

.. function:: void fmprb_gamma(fmprb_t y, const fmprb_t x, long prec)

.. function:: void fmprb_rgamma(fmprb_t y, const fmprb_t x, long prec)

.. function:: void fmprb_lgamma(fmprb_t y, const fmprb_t x, long prec)

    Sets, respectively, `y = \Gamma(x)`, `y = 1/\Gamma(x)`,
    `y = \log \Gamma(x)`. These functions are simple wrappers for the
    Stirling series code in the *gamma* module.

