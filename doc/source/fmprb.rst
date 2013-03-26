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

.. function:: void fmprb_add_error_2exp_fmpz(fmprb_t x, const fmpz_t e)

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

.. function:: long fmprb_bits(const fmprb_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of the midpoint of *x*, i.e. the minimum precision
    sufficient to represent *x* exactly. Returns 0 if the midpoint
    of *x* is a special value.


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

.. function:: void fmprb_mul_2exp_fmpz(fmprb_t y, const fmprb_t x, const fmpz_t e)

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

.. function:: void fmprb_hypot(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

    Sets *z* to `\sqrt{x^2 + y^2}`.

.. function:: void fmprb_root(fmprb_t z, const fmprb_t x, ulong k, long prec)

    Sets *z* to the *k*-th root of *x*, rounded to *prec* bits.
    As currently implemented, this function is only fast for small
    fixed *k*. For large *k* it is better to use :func:`fmprb_pow_fmpq`
    or :func:`fmprb_pow`.

.. function:: void fmprb_pow_fmpz_binexp(fmprb_t y, const fmprb_t b, const fmpz_t e, long prec)

.. function:: void fmprb_pow_fmpz(fmprb_t y, const fmprb_t b, const fmpz_t e, long prec)

.. function:: void fmprb_pow_ui(fmprb_t y, const fmprb_t b, ulong e, long prec)

.. function:: void fmprb_ui_pow_ui(fmprb_t y, ulong b, ulong e, long prec)

.. function:: void fmprb_si_pow_ui(fmprb_t y, long b, ulong e, long prec)

    Sets `y = b^e` using binary exponentiation (with an initial division
    if `e < 0`). Provided that *b* and *e*
    are small enough and the exponent is positive, the exact power can be
    computed by setting the precision to *FMPR_PREC_EXACT*.

    Note that these functions can get slow if the exponent is
    extremely large (in such cases :func:`fmprb_pow` may be superior).

.. function:: void fmprb_pow_fmpq(fmprb_t y, const fmprb_t b, const fmpq_t e, long prec)

    Sets `y = b^e`, computed as `y = (b^{1/q})^p` if the denominator of
    `e = p/q` is small, and generally as `y = \exp(e \log b)`.

    Note that this function can get slow if the exponent is
    extremely large (in such cases :func:`fmprb_pow` may be superior).

.. function:: void fmprb_pow(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

    Sets `z = x^y`, computed using binary exponentiation if `y` if
    a small exact integer, as `z = (x^{1/2})^{2y}` if `y` is a small exact
    half-integer, and generally as `z = \exp(y \log x)`.

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

.. function:: void fmprb_exp(fmprb_t z, const fmprb_t x, long prec)

    Sets `z = \exp(x)`. Error propagation is done using the following rule:
    the error is largest at `m + r`, and we have
    `\exp(m+r) - \exp(m) = \exp(m) (\exp(r)-1) \le r \exp(m+r)`.

.. function:: void fmprb_expm1(fmprb_t z, const fmprb_t x, long prec)

    Sets `z = \exp(x)-1`, computed accurately when `x \approx 0`.

Trigonometric functions
-------------------------------------------------------------------------------

.. function:: void fmprb_sin(fmprb_t s, const fmprb_t x, long prec)

.. function:: void fmprb_cos(fmprb_t c, const fmprb_t x, long prec)

.. function:: void fmprb_sin_cos(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)

    Sets `s = \sin x`, `c = \cos x`. Error propagation uses the rule
    `|\sin(m \pm r) - \sin(m)| \le \min(r,2)`.

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

    Sets `s = \sinh x`, `c = \cosh x`. If the midpoint of `x` is close
    to zero and the hyperbolic sine is to be computed,
    evaluates `(e^{2x}\pm1) / (2e^x)` via :func:`fmprb_expm1`
    to avoid loss of accuracy. Otherwise evaluates `(e^x \pm e^{-x}) / 2`.

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

    Sets *f* to the Fibonacci number `F_n`. Uses the binary squaring
    algorithm described in [Tak2000]_.
    Provided that *n* is small enough, an exact Fibonacci number can be
    computed using *FMPR_PREC_EXACT*.


Constants
-------------------------------------------------------------------------------

.. function:: void fmprb_const_pi(fmprb_t x, long prec)

    Sets *x* to `\pi`. The value is cached for repeated use.
    Uses the generic hypergeometric series code to evaluate the Chudnovsky series

    .. math ::

        \frac{1}{\pi} = 12 \sum^\infty_{k=0} \frac{(-1)^k (6k)! (13591409 + 545140134k)}{(3k)!(k!)^3 640320^{3k + 3/2}}

.. function:: void fmprb_const_sqrt_pi(fmprb_t x, long prec)

    Sets *x* to `\sqrt{\pi}`. The value is cached for repeated use.

.. function:: void fmprb_const_log2(fmprb_t s, long prec)

    Sets *x* to `\log 2`. The value is cached for repeated use.
    Uses the generic hypergeometric series code to evaluate the representation

    .. math ::

        \log 2 = \frac{3}{4} \sum_{k=0}^{\infty} \frac{(-1)^k (k!)^2}{2^k (2k+1)!}

.. function:: void fmprb_const_log10(fmprb_t x, long prec)

    Sets *x* to `\log 10`. The value is cached for repeated use.
    Uses the generic hypergeometric series code to evaluate the
    Machin-like formula
    `\log 10 = 46 \operatorname{atanh}(1/31) + 34 \operatorname{atanh}(1/49) + 20 \operatorname{atanh}(1/161)`.

.. function:: void fmprb_const_e(fmprb_t x, long prec)

    Sets *x* to Euler's number `e = \sum_{n=0}^{\infty} 1/n!`, evaluated
    using the generic hypergeometric series code.
    The value is cached for repeated use.

.. function:: void fmprb_const_euler(fmprb_t x, long prec)

    Sets *x* to Euler's constant `\gamma = \lim_{k \rightarrow \infty} (H_k - \log k)`
    where `H_k = 1 + 1/2 + \ldots + 1/k`. The value is cached for repeated use.
    Uses the Brent-McMillan formula ([BM1980]_,  [MPFR2012]_)

    .. math ::

        \gamma = \frac{S_0(2n) - K_0(2n)}{I_0(2n)} - \log(n)

    in which `n` is a free parameter and

    .. math ::

        S_0(x) = \sum_{k=0}^{\infty} \frac{H_k}{(k!)^2} \left(\frac{x}{2}\right)^{2k}, \quad
        I_0(x) = \sum_{k=0}^{\infty} \frac{1}{(k!)^2} \left(\frac{x}{2}\right)^{2k}

        2x I_0(x) K_0(x) \sim \sum_{k=0}^{\infty} \frac{[(2k)!]^3}{(k!)^4 8^{2k} x^{2k}}.

    All series are evaluated using binary splitting.
    The first two series are evaluated simultaneously, with the summation
    taken up to `\beta n` where `\beta \approx 4.9706257595442318644`
    satisfies `\beta (\log \beta - 1) = 3`. The error is easily bounded
    by twice the first discarded term and is of order `O(e^{-6n})`.
    Since both `S_0` and `I_0` are of size `O(e^{2n})`, their relative
    errors and the absolute error of their quotient is `O(e^{-8n})`.
    The third series is a divergent asymptotic expansion. With some work, it
    can be shown (to be published) that the error when `x = 2n` and
    the sum goes up to `k = 2n-1` is bounded by `8e^{-4n}`.
    Since `I_0(2n) \sim e^{2n} / n^{1/2}`, the final error is `O(e^{-8n})`.

.. function:: void fmprb_const_catalan(fmprb_t x, long prec)

    Sets *x* to Catalan's constant `C = \sum_{n=0}^{\infty} (-1)^n / (2n+1)^2`.
    The value is cached for repeated use. Uses the generic hypergeometric
    series code to evaluate the representation

    .. math ::

        C = \sum_{k=0}^{\infty} \frac{(-1)^k 4^{4 k+1}
            \left(40 k^2+56 k+19\right) [(k+1)!]^2 [(2k+2)!]^3}{(k+1)^3 (2 k+1) [(4k+4)!]^2}

.. function:: void fmprb_const_log_sqrt2pi(fmprb_t x, long prec)

    Sets *x* to `\log \sqrt{2 \pi}`. The value is cached for repeated use.


Gamma function
-------------------------------------------------------------------------------

.. function:: void fmprb_gamma(fmprb_t y, const fmprb_t x, long prec)

.. function:: void fmprb_gamma_fmpq(fmprb_t y, const fmpq_t x, long prec)

    Sets `y = \Gamma(x)`, the gamma function.

.. function:: void fmprb_rgamma(fmprb_t y, const fmprb_t x, long prec)

    Sets  `y = 1/\Gamma(x)`, avoiding division by zero at the poles
    of the gamma function.

.. function:: void fmprb_lgamma(fmprb_t y, const fmprb_t x, long prec)

    Sets `y = \log \Gamma(x)`. The complex branch structure is assumed,
    so if `x \le 0`, the result is an indeterminate interval.

.. function:: void fmprb_digamma(fmprb_t y, const fmprb_t x, long prec)

    Sets `y = \psi(x) = (\log \Gamma(x))' = \Gamma'(x) / \Gamma(x)`.

