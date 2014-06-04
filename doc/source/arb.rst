.. _arb:

**arb.h** -- real numbers represented as floating-point balls
===============================================================================

The :type:`arb_t` type is essentially identical semantically to
the :type:`fmprb_t` type, but uses an internal representation that
generally allows operation to be performed more efficiently.

Whereas the midpoint and radius of an :type:`fmprb_t` both have the
same type, the :type:`arb_t` type uses an :type:`arf_t` for the midpoint
and a :type:`mag_t` for the radius.  Code designed to manipulate the
radius of an :type:`fmprb_t` directly can be ported to the :type:`arb_t` type
by writing the radius to a temporary :type:`arf_t` variable, manipulating
that variable, and then converting back to the :type:`mag_t` radius.
Alternatively, :type:`mag_t` methods can be used directly where available.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: arb_struct

.. type:: arb_t

    An :type:`arb_struct` consists of an :type:`arf_struct` (the midpoint) and
    a :type:`mag_struct` (the radius).
    An :type:`arb_t` is defined as an array of length one of type
    :type:`arb_struct`, permitting an :type:`arb_t` to be passed by
    reference.

.. type:: arb_ptr

   Alias for ``arb_struct *``, used for vectors of numbers.

.. type:: arb_srcptr

   Alias for ``const arb_struct *``, used for vectors of numbers
   when passed as constant input to functions.

.. macro:: arb_midref(x)

    Macro returning a pointer to the midpoint of *x* as an :type:`arf_t`.

.. macro:: arb_radref(x)

    Macro returning a pointer to the radius of *x* as a :type:`mag_t`.

Memory management
-------------------------------------------------------------------------------

.. function:: void arb_init(arb_t x)

    Initializes the variable *x* for use. Its midpoint and radius are both
    set to zero.

.. function:: void arb_clear(arb_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: arb_ptr _arb_vec_init(long n)

    Returns a pointer to an array of *n* initialized :type:`arb_struct`
    entries.

.. function:: void _arb_vec_clear(arb_ptr v, long n)

    Clears an array of *n* initialized :type:`arb_struct` entries.

.. function:: void arb_swap(arb_t x, arb_t y)

    Swaps *x* and *y* efficiently.

Assignment and rounding
-------------------------------------------------------------------------------

.. function:: void arb_set_fmprb(arb_t y, const fmprb_t x)

.. function:: void arb_get_fmprb(fmprb_t y, const arb_t x)

.. function:: void arb_set(arb_t y, const arb_t x)

.. function:: void arb_set_arf(arb_t y, const arf_t x)

.. function:: void arb_set_si(arb_t y, long x)

.. function:: void arb_set_ui(arb_t y, ulong x)

.. function:: void arb_set_fmpz(arb_t y, const fmpz_t x)

    Sets *y* to the value of *x* without rounding.

.. function:: void arb_set_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t e)

    Sets *y* to `x \cdot 2^e`.

.. function:: void arb_set_round(arb_t y, const arb_t x, long prec)

.. function:: void arb_set_round_fmpz(arb_t y, const fmpz_t x, long prec)

    Sets *y* to the value of *x*, rounded to *prec* bits.

.. function:: void arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t e, long prec)

    Sets *y* to `x \cdot 2^e`, rounded to *prec* bits.

.. function:: void arb_set_fmpq(arb_t y, const fmpq_t x, long prec)

    Sets *y* to the rational number *x*, rounded to *prec* bits.

Assignment of special values
-------------------------------------------------------------------------------

.. function:: void arb_zero(arb_t x)

    Sets *x* to zero.

.. function:: void arb_one(arb_t f)

    Sets *x* to the exact integer 1.

.. function:: void arb_pos_inf(arb_t x)

    Sets *x* to positive infinity, with a zero radius.

.. function:: void arb_neg_inf(arb_t x)

    Sets *x* to negative infinity, with a zero radius.

.. function:: void arb_zero_pm_inf(arb_t x)

    Sets *x* to `[0 \pm \infty]`, representing the whole extended real line.

.. function:: void arb_indeterminate(arb_t x)

    Sets *x* to `[\operatorname{NaN} \pm \infty]`, representing
    an indeterminate result.

Input and output
-------------------------------------------------------------------------------

.. function:: void arb_print(const arb_t x)

    Prints the internal representation of *x*.

.. function:: void arb_printd(const arb_t x, long digits)

    Prints *x* in decimal. The printed value of the radius is not adjusted
    to compensate for the fact that the binary-to-decimal conversion
    of both the midpoint and the radius introduces additional error.

Random number generation
-------------------------------------------------------------------------------

.. function:: void arb_randtest(arb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random ball. The midpoint and radius will both be finite.

.. function:: void arb_randtest_exact(arb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random number with zero radius.

.. function:: void arb_randtest_precise(arb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random number with radius around `2^{-\text{prec}}`
    the magnitude of the midpoint.

.. function:: void arb_randtest_wide(arb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random number with midpoint and radius chosen independently,
    possibly giving a very large interval.

.. function:: void arb_randtest_special(arb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random interval, possibly having NaN or an infinity
    as the midpoint and possibly having an infinite radius.

.. function:: void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, long bits)

    Sets *q* to a random rational number from the interval represented by *x*.
    A denominator is chosen by multiplying the binary denominator of *x*
    by a random integer up to *bits* bits.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the exponent of the midpoint or radius is so large or small
    that representing the endpoints as exact rational numbers would
    cause overflows.

Radius and interval operations
-------------------------------------------------------------------------------

.. function:: void arb_add_error_arf(arb_t x, const arf_t err)

    Adds *err*, which is assumed to be nonnegative, to the radius of *x*.

.. function:: void arb_add_error_2exp_si(arb_t x, long e)

.. function:: void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t e)

    Adds `2^e` to the radius of *x*.

.. function:: void arb_add_error(arb_t x, const arb_t error)

    Adds the supremum of *err*, which is assumed to be nonnegative, to the
    radius of *x*.

.. function:: void arb_union(arb_t z, const arb_t x, const arb_t y, long prec)

    Sets *z* to a ball containing both *x* and *y*.

.. function:: void arb_get_abs_ubound_arf(arf_t u, const arb_t x, long prec)

    Sets *u* to the upper bound for the absolute value of *x*,
    rounded up to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void arb_get_abs_lbound_arf(arf_t u, const arb_t x, long prec)

    Sets *u* to the lower bound for the absolute value of *x*,
    rounded down to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void arb_get_mag(mag_t z, const arb_t x)

    Sets *z* to an upper bound for the absolute value of *x*. If *x* contains
    NaN, the result is positive infinity.

.. function:: void arb_get_mag_lower(mag_t z, const arb_t x)

    Sets *z* to a lower bound for the absolute value of *x*. If *x* contains
    NaN, the result is zero.

.. function:: void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x)

    Computes the exact interval represented by *x*, in the form of an integer
    interval multiplied by a power of two, i.e. `x = [a, b] \times 2^{\text{exp}}`.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the difference in magnitude between the midpoint and radius
    is so large that representing the endpoints exactly would cause overflows.

.. function:: void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, long prec)

    Sets *x* to a ball containing the interval `[a, b]`. We
    require that `a \le b`.

.. function:: long arb_rel_error_bits(const arb_t x)

    Returns the effective relative error of *x* measured in bits, defined as
    the difference between the position of the top bit in the radius
    and the top bit in the midpoint, plus one.
    The result is clamped between plus/minus *ARF_PREC_EXACT*.

.. function:: long arb_rel_accuracy_bits(const arb_t x)

    Returns the effective relative accuracy of *x* measured in bits,
    equal to the negative of the return value from :func:`arb_rel_error_bits`.

.. function:: long arb_bits(const arb_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of the midpoint of *x*, i.e. the minimum precision
    sufficient to represent *x* exactly. Returns 0 if the midpoint
    of *x* is a special value.

.. function:: void arb_trim(arb_t y, const arb_t x)

    Sets *y* to a trimmed copy of *x*: rounds *x* to a number of bits
    equal to the accuracy of *x* (as indicated by its radius),
    plus a few guard bits. The resulting ball is guaranteed to
    contain *x*, but is more economical if *x* has
    less than full accuracy.

.. function:: int arb_get_unique_fmpz(fmpz_t z, const arb_t x)

    If *x* contains a unique integer, sets *z* to that value and returns
    nonzero. Otherwise (if *x* represents no integers or more than one integer),
    returns zero.

Comparisons
-------------------------------------------------------------------------------

.. function:: int arb_is_zero(const arb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both zero.

.. function:: int arb_is_nonzero(const arb_t x)

    Returns nonzero iff zero is not contained in the interval represented
    by *x*.

.. function:: int arb_is_one(const arb_t f)

    Returns nonzero iff *x* is exactly 1.

.. function:: int arb_is_finite(const arb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both finite
    floating-point numbers, i.e. not infinities or NaN.

.. function:: int arb_is_exact(const arb_t x)

    Returns nonzero iff the radius of *x* is zero.

.. function:: int arb_is_int(const arb_t x)

    Returns nonzero iff *x* is an exact integer.

.. function:: int arb_equal(const arb_t x, const arb_t y)

    Returns nonzero iff *x* and *y* are equal as balls, i.e. have both the
    same midpoint and radius.

    Note that this is not the same thing as testing whether both
    *x* and *y* certainly represent the same real number, unless
    either *x* or *y* is exact (and neither contains NaN).
    To test whether both operands *might* represent the same mathematical
    quantity, use :func:`arb_overlaps` or :func:`arb_contains`,
    depending on the circumstance.

.. function:: int arb_is_positive(const arb_t x)

.. function:: int arb_is_nonnegative(const arb_t x)

.. function:: int arb_is_negative(const arb_t x)

.. function:: int arb_is_nonpositive(const arb_t x)

    Returns nonzero iff all points *p* in the interval represented by *x*
    satisfy, respectively, `p > 0`, `p \ge 0`, `p < 0`, `p \le 0`.
    If *x* contains NaN, returns zero.

.. function:: int arb_overlaps(const arb_t x, const arb_t y)

    Returns nonzero iff *x* and *y* have some point in common.
    If either *x* or *y* contains NaN, this function always returns nonzero
    (as a NaN could be anything, it could in particular contain any
    number that is included in the other operand).

.. function:: int arb_contains_arf(const arb_t x, const arf_t y)

.. function:: int arb_contains_fmpq(const arb_t x, const fmpq_t y)

.. function:: int arb_contains_fmpz(const arb_t x, const fmpz_t y)

.. function:: int arb_contains_si(const arb_t x, long y)

.. function:: int arb_contains_mpfr(const arb_t x, const mpfr_t y)

.. function:: int arb_contains(const arb_t x, const arb_t y)

    Returns nonzero iff the given number (or ball) *y* is contained in
    the interval represented by *x*.

    If *x* is contains NaN, this function always returns nonzero (as it
    could represent anything, and in particular could represent all
    the points included in *y*).
    If *y* contains NaN and *x* does not, it always returns zero.

.. function:: int arb_contains_zero(const arb_t x)

.. function:: int arb_contains_negative(const arb_t x)

.. function:: int arb_contains_nonpositive(const arb_t x)

.. function:: int arb_contains_positive(const arb_t x)

.. function:: int arb_contains_nonnegative(const arb_t x)

    Returns nonzero iff there is any point *p* in the interval represented
    by *x* satisfying, respectively, `p = 0`, `p < 0`, `p \le 0`, `p > 0`, `p \ge 0`.
    If *x* contains NaN, returns nonzero.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void arb_neg(arb_t y, const arb_t x)

.. function:: void arb_neg_round(arb_t y, const arb_t x, long prec)

    Sets *y* to the negation of *x*.

.. function:: void arb_abs(arb_t x, const arb_t y)

    Sets *y* to the absolute value of *x*. No attempt is made to improve the
    interval represented by *x* if it contains zero.

.. function:: void arb_add(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_add_arf(arb_t z, const arb_t x, const arf_t y, long prec)

.. function:: void arb_add_ui(arb_t z, const arb_t x, ulong y, long prec)

.. function:: void arb_add_si(arb_t z, const arb_t x, long y, long prec)

.. function:: void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    Sets `z = x + y`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t m, const fmpz_t e, long prec)

    Sets `z = x + m \cdot 2^e`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_sub(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, long prec)

.. function:: void arb_sub_ui(arb_t z, const arb_t x, ulong y, long prec)

.. function:: void arb_sub_si(arb_t z, const arb_t x, long y, long prec)

.. function:: void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    Sets `z = x - y`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_mul(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_mul_arf(arb_t z, const arb_t x, const arf_t y, long prec)

.. function:: void arb_mul_si(arb_t z, const arb_t x, long y, long prec)

.. function:: void arb_mul_ui(arb_t z, const arb_t x, ulong y, long prec)

.. function:: void arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    Sets `z = x \cdot y`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_mul_2exp_si(arb_t y, const arb_t x, long e)

.. function:: void arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)

    Sets *y* to *x* multiplied by `2^e`.

.. function:: void arb_addmul(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, long prec)

.. function:: void arb_addmul_si(arb_t z, const arb_t x, long y, long prec)

.. function:: void arb_addmul_ui(arb_t z, const arb_t x, ulong y, long prec)

.. function:: void arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    Sets `z = z + x \cdot y`, rounded to prec bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_submul(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_submul_arf(arb_t z, const arb_t x, const arf_t y, long prec)

.. function:: void arb_submul_si(arb_t z, const arb_t x, long y, long prec)

.. function:: void arb_submul_ui(arb_t z, const arb_t x, ulong y, long prec)

.. function:: void arb_submul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    Sets `z = z - x \cdot y`, rounded to prec bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_inv(arb_t y, const arb_t x, long prec)

    Sets *z* to `1 / x`.

.. function:: void arb_div(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_div_arf(arb_t z, const arb_t x, const arf_t y, long prec)

.. function:: void arb_div_si(arb_t z, const arb_t x, long y, long prec)

.. function:: void arb_div_ui(arb_t z, const arb_t x, ulong y, long prec)

.. function:: void arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

.. function:: void arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, long prec)

.. function:: void arb_ui_div(arb_t z, ulong x, const arb_t y, long prec)

    Sets `z = x / y`, rounded to *prec* bits. If *y* contains zero, *z* is
    set to `0 \pm \infty`. Otherwise, error propagation uses the rule

    .. math ::
        \left| \frac{x}{y} - \frac{x+\xi_1 a}{y+\xi_2 b} \right| =
        \left|\frac{x \xi_2 b - y \xi_1 a}{y (y+\xi_2 b)}\right| \le
        \frac{|xb|+|ya|}{|y| (|y|-b)}

    where `-1 \le \xi_1, \xi_2 \le 1`, and
    where the triangle inequality has been applied to the numerator and
    the reverse triangle inequality has been applied to the denominator.

.. function:: void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, long prec)

    Sets `z = x / (2^n - 1)`, rounded to *prec* bits.

Powers and roots
-------------------------------------------------------------------------------

.. function:: void arb_sqrt(arb_t z, const arb_t x, long prec)

.. function:: void arb_sqrt_arf(arb_t z, const arf_t x, long prec)

.. function:: void arb_sqrt_fmpz(arb_t z, const fmpz_t x, long prec)

.. function:: void arb_sqrt_ui(arb_t z, ulong x, long prec)

    Sets *z* to the square root of *x*, rounded to *prec* bits.

    If `x = m \pm x` where `m \ge r \ge 0`, the propagated error is bounded by
    `\sqrt{m} - \sqrt{m-r} = \sqrt{m} (1 - \sqrt{1 - r/m}) \le \sqrt{m} (r/m + (r/m)^2)/2`.

.. function:: void arb_sqrtpos(arb_t z, const arb_t x, long prec)

    Sets *z* to the square root of *x*, assuming that *x* represents a
    nonnegative number (i.e. discarding any negative numbers in the input
    interval), and producing an output interval not containing any
    negative numbers (unless the radius is infinite).

.. function:: void arb_hypot(arb_t z, const arb_t x, const arb_t y, long prec)

    Sets *z* to `\sqrt{x^2 + y^2}`.

.. function:: void arb_rsqrt(arb_t z, const arb_t x, long prec)

.. function:: void arb_rsqrt_ui(arb_t z, ulong x, long prec)

    Sets *z* to the reciprocal square root of *x*, rounded to *prec* bits.
    At high precision, this is faster than computing a square root.

.. function:: void arb_root(arb_t z, const arb_t x, ulong k, long prec)

    Sets *z* to the *k*-th root of *x*, rounded to *prec* bits.
    As currently implemented, this function is only fast for small *k*.
    For large *k* it is better to use :func:`arb_pow_fmpq` or :func:`arb_pow`.

.. function:: void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, long prec)

.. function:: void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, long prec)

.. function:: void arb_pow_ui(arb_t y, const arb_t b, ulong e, long prec)

.. function:: void arb_ui_pow_ui(arb_t y, ulong b, ulong e, long prec)

.. function:: void arb_si_pow_ui(arb_t y, long b, ulong e, long prec)

    Sets `y = b^e` using binary exponentiation (with an initial division
    if `e < 0`). Provided that *b* and *e*
    are small enough and the exponent is positive, the exact power can be
    computed by setting the precision to *ARF_PREC_EXACT*.

    Note that these functions can get slow if the exponent is
    extremely large (in such cases :func:`arb_pow` may be superior).

.. function:: void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, long prec)

    Sets `y = b^e`, computed as `y = (b^{1/q})^p` if the denominator of
    `e = p/q` is small, and generally as `y = \exp(e \log b)`.

    Note that this function can get slow if the exponent is
    extremely large (in such cases :func:`arb_pow` may be superior).

.. function:: void arb_pow(arb_t z, const arb_t x, const arb_t y, long prec)

    Sets `z = x^y`, computed using binary exponentiation if `y` if
    a small exact integer, as `z = (x^{1/2})^{2y}` if `y` is a small exact
    half-integer, and generally as `z = \exp(y \log x)`.

Exponentials and logarithms
-------------------------------------------------------------------------------

The following are currently simply wrappers for the corresponding *fmprb* functions.

.. function:: void arb_log(arb_t z, const arb_t x, long prec)

.. function:: void arb_log_ui(arb_t z, ulong x, long prec)

.. function:: void arb_log_fmpz(arb_t z, const fmpz_t x, long prec)

    Sets `z = \log(x)`. Error propagation is done using the following rule:
    assuming `x = m \pm r` where `m > r \ge 0`, the error is largest at
    `m - r`, and we have `\log(m) - \log(m-r) = \log(1 + r/(m-r))`.

.. function:: void arb_exp(arb_t z, const arb_t x, long prec)

    Sets `z = \exp(x)`. Error propagation is done using the following rule:
    assuming `x = m \pm r`, the error is largest at `m + r`, and we have
    `\exp(m+r) - \exp(m) = \exp(m) (\exp(r)-1) \le r \exp(m+r)`.

.. function:: void arb_expm1(arb_t z, const arb_t x, long prec)

    Sets `z = \exp(x)-1`, computed accurately when `x \approx 0`.

Trigonometric functions
-------------------------------------------------------------------------------

The following are currently simply wrappers for the corresponding *fmprb* functions.

.. function:: void arb_sin(arb_t s, const arb_t x, long prec)

.. function:: void arb_cos(arb_t c, const arb_t x, long prec)

.. function:: void arb_sin_cos(arb_t s, arb_t c, const arb_t x, long prec)

    Sets `s = \sin x`, `c = \cos x`. Error propagation uses the rule
    `|\sin(m \pm r) - \sin(m)| \le \min(r,2)`.

.. function:: void arb_sin_pi(arb_t s, const arb_t x, long prec)

.. function:: void arb_cos_pi(arb_t c, const arb_t x, long prec)

.. function:: void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, long prec)

    Sets `s = \sin \pi x`, `c = \cos \pi x`.

.. function:: void arb_tan(arb_t y, const arb_t x, long prec)

    Sets `y = \tan x = (\sin x) / (\cos y)`.

.. function:: void arb_cot(arb_t y, const arb_t x, long prec)

    Sets `y = \cot x = (\cos x) / (\sin y)`.

.. function:: void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, long prec)

.. function:: void arb_sin_pi_fmpq(arb_t s, const fmpq_t x, long prec)

.. function:: void arb_cos_pi_fmpq(arb_t c, const fmpq_t x, long prec)

    Sets `s = \sin \pi x`, `c = \cos \pi x` where `x` is a rational
    number (whose numerator and denominator are assumed to be reduced).
    We first use trigonometric symmetries to reduce the argument to the
    octant `[0, 1/4]`. Then we either multiply by a numerical approximation
    of `\pi` and evaluate the trigonometric function the usual way,
    or we use algebraic methods, depending on which is estimated to be faster.
    Since the argument has been reduced to the first octant, the
    first of these two methods gives full accuracy even if the original
    argument is close to some root other the origin.

.. function:: void arb_tan_pi(arb_t y, const arb_t x, long prec)

    Sets `y = \tan \pi x`.

.. function:: void arb_cot_pi(arb_t y, const arb_t x, long prec)

    Sets `y = \cot \pi x`.

Inverse trigonometric functions
-------------------------------------------------------------------------------

The following are currently simply wrappers for the corresponding *fmprb* functions.

.. function:: void arb_atan(arb_t z, const arb_t x, long prec)

    Sets `z = \tan^{-1} x`. Letting `d = \max(0, |m| - r)`,
    the propagated error is bounded by `r / (1 + d^2)`
    (this could be tightened).

.. function:: void arb_atan2(arb_t z, const arb_t b, const arb_t a, long prec)

    Sets *r* to an the argument (phase) of the complex number
    `a + bi`, with the branch cut discontinuity on `(-\infty,0]`.
    We define `\operatorname{atan2}(0,0) = 0`, and for `a < 0`,
    `\operatorname{atan2}(0,a) = \pi`.

.. function:: void arb_asin(arb_t z, const arb_t x, long prec)

    Sets `z = \sin^{-1} x = \tan^{-1}(x / \sqrt{1-x^2})`.
    If `x` is not contained in the domain `[-1,1]`, the result is an
    indeterminate interval.

.. function:: void arb_acos(arb_t z, const arb_t x, long prec)

    Sets `z = \cos^{-1} x = \pi/2 - \sin^{-1} x`.
    If `x` is not contained in the domain `[-1,1]`, the result is an
    indeterminate interval.

Hyperbolic functions
-------------------------------------------------------------------------------

.. function:: void arb_sinh(arb_t s, const arb_t x, long prec)

.. function:: void arb_cosh(arb_t c, const arb_t x, long prec)

.. function:: void arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, long prec)

    Sets `s = \sinh x`, `c = \cosh x`. If the midpoint of `x` is close
    to zero and the hyperbolic sine is to be computed,
    evaluates `(e^{2x}\pm1) / (2e^x)` via :func:`arb_expm1`
    to avoid loss of accuracy. Otherwise evaluates `(e^x \pm e^{-x}) / 2`.

.. function:: void arb_tanh(arb_t y, const arb_t x, long prec)

    Sets `y = \tanh x = (\sinh x) / (\cosh x)`, evaluated
    via :func:`arb_expm1` as `\tanh x = (e^{2x} - 1) / (e^{2x} + 1)` if
    the midpoint of `x` is negative and as
    `\tanh x = (1 - e^{-2x}) / (1 + e^{-2x})` otherwise.

.. function:: void arb_coth(arb_t y, const arb_t x, long prec)

    Sets `y = \coth x = (\cosh x) / (\sinh x)`, evaluated using
    the same strategy as :func:`arb_tanh`.

Constants
-------------------------------------------------------------------------------

The following are currently simply wrappers for the corresponding *fmprb* functions.

.. function:: void arb_const_pi(arb_t z, long prec)

.. function:: void arb_const_sqrt_pi(arb_t z, long prec)

.. function:: void arb_const_log_sqrt2pi(arb_t z, long prec)

.. function:: void arb_const_log2(arb_t z, long prec)

.. function:: void arb_const_log10(arb_t z, long prec)

.. function:: void arb_const_euler(arb_t z, long prec)

.. function:: void arb_const_catalan(arb_t z, long prec)

.. function:: void arb_const_e(arb_t z, long prec)

.. function:: void arb_const_khinchin(arb_t z, long prec)

.. function:: void arb_const_glaisher(arb_t z, long prec)

Rising factorials
-------------------------------------------------------------------------------

.. function:: void arb_rising_ui_bs(arb_t z, const arb_t x, ulong n, long prec)

.. function:: void arb_rising_ui_rs(arb_t z, const arb_t x, ulong n, ulong step, long prec)

.. function:: void arb_rising_ui_rec(arb_t z, const arb_t x, ulong n, long prec)

.. function:: void arb_rising_ui(arb_t z, const arb_t x, ulong n, long prec)

    Computes the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)`.

    The *bs* version uses binary splitting. The *rs* version uses rectangular
    splitting. The *rec* version uses either *bs* or *rs* depending
    on the input.
    The default version is currently identical to the *rec* version.
    In a future version, it will use the gamma function or asymptotic
    series when this is more efficient.

    The *rs* version takes an optional *step* parameter for tuning
    purposes (to use the default step length, pass zero).

.. function:: void arb_rising_fmpq_ui(arb_t z, const fmpq_t x, ulong n, long prec)

    Computes the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)` using
    binary splitting. If the denominator or numerator of *x* is large
    compared to *prec*, it is more efficient to convert *x* to an approximation
    and use :func:`arb_rising_ui`.

Special functions
-------------------------------------------------------------------------------

The following are currently simply wrappers for the corresponding *fmprb* functions.

.. function:: void arb_fac_ui(arb_t z, ulong n, long prec)

.. function:: void arb_bin_ui(arb_t z, const arb_t n, ulong k, long prec)

.. function:: void arb_bin_uiui(arb_t z, ulong n, ulong k, long prec)

.. function:: void arb_fib_fmpz(arb_t z, const fmpz_t n, long prec)

.. function:: void arb_fib_ui(arb_t z, ulong n, long prec)

.. function:: void arb_agm(arb_t z, const arb_t x, const arb_t y, long prec)

.. function:: void arb_lgamma(arb_t z, const arb_t x, long prec)

.. function:: void arb_rgamma(arb_t z, const arb_t x, long prec)

.. function:: void arb_gamma(arb_t z, const arb_t x, long prec)

.. function:: void arb_gamma_fmpq(arb_t z, const fmpq_t x, long prec)

.. function:: void arb_gamma_fmpz(arb_t z, const fmpz_t x, long prec)

.. function:: void arb_digamma(arb_t y, const arb_t x, long prec)

.. function:: void arb_zeta(arb_t z, const arb_t s, long prec)

.. function:: void arb_zeta_ui(arb_t z, ulong n, long prec)

.. function:: void arb_bernoulli_ui(arb_t z, ulong n, long prec)

