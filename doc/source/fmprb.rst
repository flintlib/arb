.. _fmprb:

**fmprb.h** -- real numbers represented as floating-point balls
===============================================================================

An :type:`fmprb_t` represents a ball over the real numbers,
that is, an interval `[m \pm r] \equiv [m-r, m+r]` where the midpoint `m` and the
radius `r` are (extended) real numbers and `r` is nonnegative (possibly infinite).
The result of an (approximate) operation done on *fmprb_t* variables
is a ball which contains the result of the (mathematically exact) operation
applied to any choice of points in the input balls.
In general, the output ball is not the smallest possible.

The precision (*prec*) parameter passed to each function roughly indicates the
precision to which calculations on the midpoint are carried out (operations on
the radius are always done using a fixed, small precision.) In other words, the
*prec* parameter just specifies the internal working precision used to compute
the midpoint. For example, if adding two balls `[x \pm r]`, `[y \pm s]`, the
result is
`[\operatorname{round}(x + y, \operatorname{prec}) \pm (r + s + \varepsilon)]`
where `\varepsilon` is a bound for the error from
`\operatorname{round}(x + y, \operatorname{prec})`.

If the precision is much higher than the accuracy of the inputs, all that
happens is that `\varepsilon` is much smaller than `r + s`, so the output
radius is approximately `r + s`.

If the precision is much lower than the accuracy of the inputs, all that
happens is that `\varepsilon` is much larger than `r + s`, so the output radius
is approximately `\varepsilon`.

The only restrictions are that the precision has to be at least 2 bits, and not
so large that overflow occurs (not usually an issue, at least on 64-bit
systems).

The only strict guarantee for all functions is that the output ball contains
the mathematically exact result.

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

The following balls with an infinite or NaN component are permitted,
and may be returned as output from functions.

* The ball `[+\infty \pm c]`, where `c` is finite, represents the point at positive infinity. Such a ball can always be replaced by `[+\infty \pm 0]` while preserving mathematical correctness (this is currently not done automatically by the library).
* The ball `[-\infty \pm c]`, where `c` is finite, represents the point at negative infinity. Such a ball can always be replaced by `[-\infty \pm 0]` while preserving mathematical correctness (this is currently not done automatically by the library).
* The ball `[c \pm \infty]`, where `c` is finite or infinite, represents the whole extended real line `[-\infty,+\infty]`. Such a ball can always be replaced by `[0 \pm \infty]` while preserving mathematical correctness (this is currently not done automatically by the library). Note that there is no way to represent a half-infinite interval such as `[0,\infty]`.
* The ball `[\operatorname{NaN} \pm c]`, where `c` is finite or infinite, represents an indeterminate value (the value could be any extended real number, or it could represent a function being evaluated outside its domain of definition, for example where the result would be complex). Such an indeterminate ball can always be replaced by `[\operatorname{NaN} \pm \infty]` while preserving mathematical correctness (this is currently not done automatically by the library).

The radius of a ball is not allowed to be negative or NaN.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmprb_struct

.. type:: fmprb_t

    An *fmprb_struct* consists of a pair of *fmpr_struct*:s.
    An *fmprb_t* is defined as an array of length one of type
    *fmprb_struct*, permitting an *fmprb_t* to be passed by
    reference.

.. type:: fmprb_ptr

   Alias for ``fmprb_struct *``, used for vectors of numbers.

.. type:: fmprb_srcptr

   Alias for ``const fmprb_struct *``, used for vectors of numbers
   when passed as constant input to functions.

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

.. function:: fmprb_ptr _fmprb_vec_init(long n)

    Returns a pointer to an array of *n* initialized *fmprb_struct*:s.

.. function:: void _fmprb_vec_clear(fmprb_ptr v, long n)

    Clears an array of *n* initialized *fmprb_struct*:s.


Assignment and rounding
-------------------------------------------------------------------------------

.. function:: void fmprb_set(fmprb_t y, const fmprb_t x)

    Sets *y* to a copy of *x*.

.. function:: void fmprb_set_round(fmprb_t y, const fmprb_t x, long prec)

    Sets *y* to a copy of *x*, rounded to *prec* bits.

.. function:: void fmprb_set_fmpr(fmprb_t y, const fmpr_t x)

.. function:: void fmprb_set_si(fmprb_t y, long x)

.. function:: void fmprb_set_ui(fmprb_t y, ulong x)

.. function:: void fmprb_set_fmpz(fmprb_t y, const fmpz_t x)

    Sets *y* exactly to *x*.

.. function:: void fmprb_set_fmpq(fmprb_t y, const fmpq_t x, long prec)

    Sets *y* to the rational number *x*, rounded to *prec* bits.

.. function:: void fmprb_set_fmpz_2exp(fmprb_t x, const fmpz_t y, const fmpz_t exp)

    Sets *x* to *y* multiplied by 2 raised to the power *exp*.

.. function:: void fmprb_set_round_fmpz_2exp(fmprb_t y, const fmpz_t x, const fmpz_t exp, long prec)

    Sets *x* to *y* multiplied by 2 raised to the power *exp*, rounding
    the result to *prec* bits.


Assignment of special values
-------------------------------------------------------------------------------

.. function:: void fmprb_zero(fmprb_t x)

    Sets *x* to zero.

.. function:: void fmprb_one(fmprb_t x)

    Sets *x* to the exact integer 1.

.. function:: void fmprb_pos_inf(fmprb_t x)

    Sets *x* to positive infinity, with a zero radius.

.. function:: void fmprb_neg_inf(fmprb_t x)

    Sets *x* to negative infinity, with a zero radius.

.. function:: void fmprb_zero_pm_inf(fmprb_t x)

    Sets *x* to `[0 \pm \infty]`, representing the whole extended real line.

.. function:: void fmprb_indeterminate(fmprb_t x)

    Sets *x* to `[\operatorname{NaN} \pm \infty]`, representing
    an indeterminate result.


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

.. function:: void fmprb_randtest_special(fmprb_t x, flint_rand_t state, long prec, long mag_bits)

    Generates a random interval, possibly having NaN or an infinity
    as the midpoint and possibly having an infinite radius.

.. function:: void fmprb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const fmprb_t x, long bits)

    Sets *q* to a random rational number from the interval represented by *x*.
    A denominator is chosen by multiplying the binary denominator of *x*
    by a random integer up to *bits* bits.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the exponent of the midpoint or radius is so large or small
    that representing the endpoints as exact rational numbers would
    cause overflows.


Radius and interval operations
-------------------------------------------------------------------------------

.. function:: void fmprb_add_error_fmpr(fmprb_t x, const fmpr_t err)

    Adds *err*, which is assumed to be nonnegative, to the radius of *x*.

.. function:: void fmprb_add_error_2exp_si(fmprb_t x, long e)

.. function:: void fmprb_add_error_2exp_fmpz(fmprb_t x, const fmpz_t e)

    Adds `2^e` to the radius of *x*.

.. function:: void fmprb_add_error(fmprb_t x, const fmprb_t err)

    Adds the supremum of *err*, which is assumed to be nonnegative, to the
    radius of *x*.

.. function:: void fmprb_union(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

    Sets *z* to a ball containing both *x* and *y*.

.. function:: void fmprb_get_abs_ubound_fmpr(fmpr_t u, const fmprb_t x, long prec)

    Sets *u* to the upper bound of the absolute value of *x*,
    rounded up to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void fmprb_get_abs_lbound_fmpr(fmpr_t u, const fmprb_t x, long prec)

    Sets *u* to the lower bound of the absolute value of *x*,
    rounded down to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void fmprb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const fmprb_t x)

    Computes the exact interval represented by *x*, in the form of an integer
    interval multiplied by a power of two, i.e. `x = [a, b] \times 2^{\mathrm{exp}}`.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the difference in magnitude between the midpoint and radius
    is so large that representing the endpoints exactly would cause overflows.

.. function:: void fmprb_set_interval_fmpr(fmprb_t x, const fmpr_t a, const fmpr_t b, long prec)

    Sets *x* to a ball containing the interval `[a, b]`. We
    require that `a \le b`.

.. function:: long fmprb_rel_error_bits(const fmprb_t x)

    Returns the effective relative error of *x* measured in bits, defined as
    the difference between the position of the top bit in the radius
    and the top bit in the midpoint, plus one.
    The result is clamped between plus/minus *FMPR_PREC_EXACT*.

.. function:: long fmprb_rel_accuracy_bits(const fmprb_t x)

    Returns the effective relative accuracy of *x* measured in bits,
    equal to the negative of the return value from *fmprb_rel_error_bits*.

.. function:: long fmprb_bits(const fmprb_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of the midpoint of *x*, i.e. the minimum precision
    sufficient to represent *x* exactly. Returns 0 if the midpoint
    of *x* is a special value.

.. function:: void fmprb_trim(fmprb_t y, const fmprb_t x)

    Sets *y* to a trimmed copy of *x*: rounds *x* to a number of bits
    equal to the accuracy of *x* (as indicated by its radius),
    plus a few guard bits. The resulting ball is guaranteed to
    contain *x*, but is more economical if *x* has
    less than full accuracy.

.. function:: int fmprb_get_unique_fmpz(fmpz_t z, const fmprb_t x)

    If *x* contains a unique integer, sets *z* to that value and returns
    nonzero. Otherwise (if *x* represents no integers or more than one integer),
    returns zero.


Comparisons
-------------------------------------------------------------------------------

.. function:: int fmprb_is_zero(const fmprb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both zero.

.. function:: int fmprb_is_nonzero(const fmprb_t x)

    Returns nonzero iff zero is not contained in the interval represented
    by *x*.

.. function:: int fmprb_is_one(const fmprb_t x)

    Returns nonzero iff *x* is exactly 1.

.. function:: int fmprb_is_finite(fmprb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both finite
    floating-point numbers, i.e. not infinities or NaN.

.. function:: int fmprb_is_exact(const fmprb_t x)

    Returns nonzero iff the radius of *x* is zero.

.. function:: int fmprb_is_int(const fmprb_t x)

    Returns nonzero iff *x* is an exact integer.

.. function:: int fmprb_equal(const fmprb_t x, const fmprb_t y)

    Returns nonzero iff *x* and *y* are equal as balls, i.e. have both the
    same midpoint and radius.

    Note that this is not the same thing as testing whether both
    *x* and *y* certainly represent the same real number, unless
    either *x* or *y* is exact (and neither contains NaN).
    To test whether both operands *might* represent the same mathematical
    quantity, use :func:`fmprb_overlaps` or :func:`fmprb_contains`,
    depending on the circumstance.

.. function:: int fmprb_is_positive(const fmprb_t x)

.. function:: int fmprb_is_nonnegative(const fmprb_t x)

.. function:: int fmprb_is_negative(const fmprb_t x)

.. function:: int fmprb_is_nonpositive(const fmprb_t x)

    Returns nonzero iff all points *p* in the interval represented by *x*
    satisfy, respectively, `p > 0`, `p \ge 0`, `p < 0`, `p \le 0`.
    If *x* contains NaN, returns zero.

.. function:: int fmprb_overlaps(const fmprb_t x, const fmprb_t y)

    Returns nonzero iff *x* and *y* have some point in common.
    If either *x* or *y* contains NaN, this function always returns nonzero
    (as a NaN could be anything, it could in particular contain any
    number that is included in the other operand).

.. function:: int fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y)

.. function:: int fmprb_contains_fmpq(const fmprb_t x, const fmpq_t y)

.. function:: int fmprb_contains_fmpz(const fmprb_t x, const fmpz_t y)

.. function:: int fmprb_contains_si(const fmprb_t x, long y)

.. function:: int fmprb_contains_mpfr(const fmprb_t x, const mpfr_t y)

.. function:: int fmprb_contains_zero(const fmprb_t x)

.. function:: int fmprb_contains(const fmprb_t x, const fmprb_t y)

    Returns nonzero iff the given number (or ball) *y* is contained in
    the interval represented by *x*.

    If *x* is contains NaN, this function always returns nonzero (as it
    could represent anything, and in particular could represent all
    the points included in *y*).
    If *y* contains NaN and *x* does not, it always returns zero.

.. function:: int fmprb_contains_negative(const fmprb_t x)

.. function:: int fmprb_contains_nonpositive(const fmprb_t x)

.. function:: int fmprb_contains_positive(const fmprb_t x)

.. function:: int fmprb_contains_nonnegative(const fmprb_t x)

    Returns nonzero iff there is any point *p* in the interval represented
    by *x* satisfying, respectively, `p < 0`, `p \le 0`, `p > 0`, `p \ge 0`.
    If *x* contains NaN, returns nonzero.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmprb_neg(fmprb_t y, const fmprb_t x)

    Sets *y* to the negation of *x*.

.. function:: void fmprb_abs(fmprb_t y, const fmprb_t x)

    Sets *y* to the absolute value of *x*. No attempt is made to improve the
    interval represented by *x* if it contains zero.

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

.. function:: void fmprb_inv(fmprb_t z, const fmprb_t x, long prec)

    Sets *z* to the multiplicative inverse of *x*.

.. function:: void fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

.. function:: void fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)

.. function:: void fmprb_div_si(fmprb_t z, const fmprb_t x, long y, long prec)

.. function:: void fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)

.. function:: void fmprb_div_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)

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

.. function:: void fmprb_rsqrt(fmprb_t z, const fmprb_t x, long prec)

.. function:: void fmprb_rsqrt_ui(fmprb_t z, ulong x, long prec)

    Sets *z* to the reciprocal square root of *x*, rounded to *prec* bits.
    At high precision, this is faster than computing a square root.

.. function:: void fmprb_root(fmprb_t z, const fmprb_t x, ulong k, long prec)

    Sets *z* to the *k*-th root of *x*, rounded to *prec* bits.
    As currently implemented, this function is only fast for small
    fixed *k*. For large *k* it is better to use :func:`fmprb_pow_fmpq`
    or :func:`fmprb_pow`.

.. function:: void fmprb_agm(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)

    Sets *z* to the arithmetic-geometric mean of *x* and *y*.

