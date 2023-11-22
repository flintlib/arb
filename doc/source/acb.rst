.. _acb:

**acb.h** -- complex numbers
===============================================================================

An :type:`acb_t` represents a complex number with
error bounds. An :type:`acb_t` consists of a pair of real number
balls of type :type:`arb_struct`, representing the real and
imaginary part with separate error bounds.

An :type:`acb_t` thus represents a rectangle
`[m_1-r_1, m_1+r_1] + [m_2-r_2, m_2+r_2] i` in the complex plane.
This is used instead of a disk or square representation
(consisting of a complex floating-point midpoint with a single radius),
since it allows implementing many operations more conveniently by splitting
into ball operations on the real and imaginary parts.
It also allows tracking when complex numbers have an exact (for example
exactly zero) real part and an inexact imaginary part, or vice versa.

The interface for the :type:`acb_t` type is slightly less developed
than that for the :type:`arb_t` type. In many cases, the user can
easily perform missing operations by directly manipulating the real and
imaginary parts.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: acb_struct

.. type:: acb_t

    An *acb_struct* consists of a pair of *arb_struct*:s.
    An *acb_t* is defined as an array of length one of type
    *acb_struct*, permitting an *acb_t* to be passed by
    reference.

.. type:: acb_ptr

   Alias for ``acb_struct *``, used for vectors of numbers.

.. type:: acb_srcptr

   Alias for ``const acb_struct *``, used for vectors of numbers
   when passed as constant input to functions.

.. macro:: acb_realref(x)

    Macro returning a pointer to the real part of *x* as an *arb_t*.

.. macro:: acb_imagref(x)

    Macro returning a pointer to the imaginary part of *x* as an *arb_t*.

Memory management
-------------------------------------------------------------------------------

.. function:: void acb_init(acb_t x)

    Initializes the variable *x* for use, and sets its value to zero.

.. function:: void acb_clear(acb_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: acb_ptr _acb_vec_init(slong n)

    Returns a pointer to an array of *n* initialized *acb_struct*:s.

.. function:: void _acb_vec_clear(acb_ptr v, slong n)

    Clears an array of *n* initialized *acb_struct*:s.

.. function:: slong acb_allocated_bytes(const acb_t x)

    Returns the total number of bytes heap-allocated internally by this object.
    The count excludes the size of the structure itself. Add
    ``sizeof(acb_struct)`` to get the size of the object as a whole.

.. function:: slong _acb_vec_allocated_bytes(acb_srcptr vec, slong len)

    Returns the total number of bytes allocated for this vector, i.e. the
    space taken up by the vector itself plus the sum of the internal heap
    allocation sizes for all its member elements.

.. function:: double _acb_vec_estimate_allocated_bytes(slong len, slong prec)

    Estimates the number of bytes that need to be allocated for a vector of
    *len* elements with *prec* bits of precision, including the space for
    internal limb data.
    See comments for :func:`_arb_vec_estimate_allocated_bytes`.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void acb_zero(acb_t z)

.. function:: void acb_one(acb_t z)

.. function:: void acb_onei(acb_t z)

    Sets *z* respectively to 0, 1, `i = \sqrt{-1}`.

.. function:: void acb_set(acb_t z, const acb_t x)

.. function:: void acb_set_ui(acb_t z, ulong x)

.. function:: void acb_set_si(acb_t z, slong x)

.. function:: void acb_set_d(acb_t z, double x)

.. function:: void acb_set_fmpz(acb_t z, const fmpz_t x)

.. function:: void acb_set_arb(acb_t z, const arb_t c)

    Sets *z* to the value of *x*.

.. function:: void acb_set_si_si(acb_t z, slong x, slong y)

.. function:: void acb_set_d_d(acb_t z, double x, double y)

.. function:: void acb_set_fmpz_fmpz(acb_t z, const fmpz_t x, const fmpz_t y)

.. function:: void acb_set_arb_arb(acb_t z, const arb_t x, const arb_t y)

    Sets the real and imaginary part of *z* to the values *x* and *y* respectively

.. function:: void acb_set_fmpq(acb_t z, const fmpq_t x, slong prec)

.. function:: void acb_set_round(acb_t z, const acb_t x, slong prec)

.. function:: void acb_set_round_fmpz(acb_t z, const fmpz_t x, slong prec)

.. function:: void acb_set_round_arb(acb_t z, const arb_t x, slong prec)

    Sets *z* to *x*, rounded to *prec* bits.

.. function:: void acb_swap(acb_t z, acb_t x)

    Swaps *z* and *x* efficiently.

.. function:: void acb_add_error_arf(acb_t x, const arf_t err)

.. function:: void acb_add_error_mag(acb_t x, const mag_t err)

.. function:: void acb_add_error_arb(acb_t x, const arb_t err)

    Adds *err* to the error bounds of both the real and imaginary
    parts of *x*, modifying *x* in-place.

.. function:: void acb_get_mid(acb_t m, const acb_t x)

    Sets *m* to the midpoint of *x*.

Input and output
-------------------------------------------------------------------------------

The *acb_print...* functions print to standard output, while
*acb_fprint...* functions print to the stream *file*.

.. function:: void acb_print(const acb_t x)

.. function:: void acb_fprint(FILE * file, const acb_t x)

    Prints the internal representation of *x*.

.. function:: void acb_printd(const acb_t x, slong digits)

.. function:: void acb_fprintd(FILE * file, const acb_t x, slong digits)

    Prints *x* in decimal. The printed value of the radius is not adjusted
    to compensate for the fact that the binary-to-decimal conversion
    of both the midpoint and the radius introduces additional error.

.. function:: void acb_printn(const acb_t x, slong digits, ulong flags)

.. function:: void acb_fprintn(FILE * file, const acb_t x, slong digits, ulong flags)

    Prints a nice decimal representation of *x*, using the format of
    :func:`arb_get_str` (or the corresponding :func:`arb_printn`) for the 
    real and imaginary parts.

    By default, the output shows the midpoint of both the real and imaginary
    parts with a guaranteed error of at most one unit in the last decimal
    place. In addition, explicit error bounds are printed so that the displayed
    decimal interval is guaranteed to enclose *x*.

    Any flags understood by :func:`arb_get_str` can be passed via *flags*
    to control the format of the real and imaginary parts.

Random number generation
-------------------------------------------------------------------------------

.. function:: void acb_randtest(acb_t z, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random complex number by generating separate random
    real and imaginary parts.

.. function:: void acb_randtest_special(acb_t z, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random complex number by generating separate random
    real and imaginary parts. Also generates NaNs and infinities.

.. function:: void acb_randtest_precise(acb_t z, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random complex number with precise real and imaginary parts.

.. function:: void acb_randtest_param(acb_t z, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random complex number, with very high probability of
    generating integers and half-integers.

Precision and comparisons
-------------------------------------------------------------------------------

.. function:: int acb_is_zero(const acb_t z)

    Returns nonzero iff *z* is zero.

.. function:: int acb_is_one(const acb_t z)

    Returns nonzero iff *z* is exactly 1.

.. function:: int acb_is_finite(const acb_t z)

    Returns nonzero iff *z* certainly is finite.

.. function:: int acb_is_exact(const acb_t z)

    Returns nonzero iff *z* is exact.

.. function:: int acb_is_int(const acb_t z)

    Returns nonzero iff *z* is an exact integer.

.. function:: int acb_is_int_2exp_si(const acb_t x, slong e)

    Returns nonzero iff *z* exactly equals `n 2^e` for some integer *n*.

.. function:: int acb_equal(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* are identical as sets, i.e.
    if the real and imaginary parts are equal as balls.

    Note that this is not the same thing as testing whether both
    *x* and *y* certainly represent the same complex number, unless
    either *x* or *y* is exact (and neither contains NaN).
    To test whether both operands *might* represent the same mathematical
    quantity, use :func:`acb_overlaps` or :func:`acb_contains`,
    depending on the circumstance.

.. function:: int acb_equal_si(const acb_t x, slong y)

    Returns nonzero iff *x* is equal to the integer *y*.

.. function:: int acb_eq(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* are certainly equal, as determined
    by testing that :func:`arb_eq` holds for both the real and imaginary
    parts.

.. function:: int acb_ne(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* are certainly not equal, as determined
    by testing that :func:`arb_ne` holds for either the real or imaginary parts.

.. function:: int acb_overlaps(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* have some point in common.

.. function:: void acb_union(acb_t z, const acb_t x, const acb_t y, slong prec)

    Sets *z* to a complex interval containing both *x* and *y*.

.. function:: void acb_get_abs_ubound_arf(arf_t u, const acb_t z, slong prec)

    Sets *u* to an upper bound for the absolute value of *z*, computed
    using a working precision of *prec* bits.

.. function:: void acb_get_abs_lbound_arf(arf_t u, const acb_t z, slong prec)

    Sets *u* to a lower bound for the absolute value of *z*, computed
    using a working precision of *prec* bits.

.. function:: void acb_get_rad_ubound_arf(arf_t u, const acb_t z, slong prec)

    Sets *u* to an upper bound for the error radius of *z* (the value
    is currently not computed tightly).

.. function:: void acb_get_mag(mag_t u, const acb_t x)

    Sets *u* to an upper bound for the absolute value of *x*.

.. function:: void acb_get_mag_lower(mag_t u, const acb_t x)

    Sets *u* to a lower bound for the absolute value of *x*.

.. function:: int acb_contains_fmpq(const acb_t x, const fmpq_t y)

.. function:: int acb_contains_fmpz(const acb_t x, const fmpz_t y)

.. function:: int acb_contains(const acb_t x, const acb_t y)

    Returns nonzero iff *y* is contained in *x*.

.. function:: int acb_contains_zero(const acb_t x)

    Returns nonzero iff zero is contained in *x*.

.. function:: int acb_contains_int(const acb_t x)

    Returns nonzero iff the complex interval represented by *x* contains
    an integer.

.. function:: int acb_contains_interior(const acb_t x, const acb_t y)

    Tests if *y* is contained in the interior of *x*.
    This predicate always evaluates to false if *x* and *y* are both
    real-valued, since an imaginary part of 0 is not considered contained in
    the interior of the point interval 0. More generally, the same
    problem occurs for intervals with an exact real or imaginary part.
    Such intervals must be handled specially by the user where a different
    interpretation is intended.

.. function:: slong acb_rel_error_bits(const acb_t x)

    Returns the effective relative error of *x* measured in bits.
    This is computed as if calling :func:`arb_rel_error_bits` on the
    real ball whose midpoint is the larger out of the real and imaginary
    midpoints of *x*, and whose radius is the larger out of the real
    and imaginary radiuses of *x*.

.. function:: slong acb_rel_accuracy_bits(const acb_t x)

    Returns the effective relative accuracy of *x* measured in bits,
    equal to the negative of the return value from :func:`acb_rel_error_bits`.

.. function:: slong acb_rel_one_accuracy_bits(const acb_t x)

    Given a ball with midpoint *m* and radius *r*, returns an approximation of
    the relative accuracy of `[\max(1,|m|) \pm r]` measured in bits.

.. function:: slong acb_bits(const acb_t x)

    Returns the maximum of *arb_bits* applied to the real
    and imaginary parts of *x*, i.e. the minimum precision sufficient
    to represent *x* exactly.

.. function:: void acb_indeterminate(acb_t x)

    Sets *x* to
    `[\operatorname{NaN} \pm \infty] + [\operatorname{NaN} \pm \infty]i`,
    representing an indeterminate result.

.. function:: void acb_trim(acb_t y, const acb_t x)

    Sets *y* to a a copy of *x* with both the real and imaginary
    parts trimmed (see :func:`arb_trim`).

.. function:: int acb_is_real(const acb_t x)

    Returns nonzero iff the imaginary part of *x* is zero.
    It does not test whether the real part of *x* also is finite.

.. function:: int acb_get_unique_fmpz(fmpz_t z, const acb_t x)

    If *x* contains a unique integer, sets *z* to that value and returns
    nonzero. Otherwise (if *x* represents no integers or more than one integer),
    returns zero.

Complex parts
-------------------------------------------------------------------------------

.. function:: void acb_get_real(arb_t re, const acb_t z)

    Sets *re* to the real part of *z*.

.. function:: void acb_get_imag(arb_t im, const acb_t z)

    Sets *im* to the imaginary part of *z*.

.. function:: void acb_arg(arb_t r, const acb_t z, slong prec)

    Sets *r* to a real interval containing the complex argument (phase) of *z*.
    We define the complex argument have a discontinuity on `(-\infty,0]`, with
    the special value `\operatorname{arg}(0) = 0`, and
    `\operatorname{arg}(a+0i) = \pi` for `a < 0`. Equivalently, if
    `z = a+bi`, the argument is given by `\operatorname{atan2}(b,a)`
    (see :func:`arb_atan2`).

.. function:: void acb_abs(arb_t r, const acb_t z, slong prec)

    Sets *r* to the absolute value of *z*.

.. function:: void acb_sgn(acb_t r, const acb_t z, slong prec)

    Sets *r* to the complex sign of *z*, defined as 0 if *z* is exactly zero
    and the projection onto the unit circle `z / |z| = \exp(i \arg(z))` otherwise.

.. function:: void acb_csgn(arb_t r, const acb_t z)

    Sets *r* to the extension of the real sign function taking the
    value 1 for *z* strictly in the right half plane, -1 for *z* strictly
    in the left half plane, and the sign of the imaginary part when *z* is on
    the imaginary axis. Equivalently, `\operatorname{csgn}(z) = z / \sqrt{z^2}`
    except that the value is 0 when *z* is exactly zero.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void acb_neg(acb_t z, const acb_t x)

.. function:: void acb_neg_round(acb_t z, const acb_t x, slong prec)

    Sets *z* to the negation of *x*.

.. function:: void acb_conj(acb_t z, const acb_t x)

    Sets *z* to the complex conjugate of *x*.

.. function:: void acb_add_ui(acb_t z, const acb_t x, ulong y, slong prec)

.. function:: void acb_add_si(acb_t z, const acb_t x, slong y, slong prec)

.. function:: void acb_add_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)

.. function:: void acb_add_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

.. function:: void acb_add(acb_t z, const acb_t x, const acb_t y, slong prec)

    Sets *z* to the sum of *x* and *y*.

.. function:: void acb_sub_ui(acb_t z, const acb_t x, ulong y, slong prec)

.. function:: void acb_sub_si(acb_t z, const acb_t x, slong y, slong prec)

.. function:: void acb_sub_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)

.. function:: void acb_sub_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

.. function:: void acb_sub(acb_t z, const acb_t x, const acb_t y, slong prec)

    Sets *z* to the difference of *x* and *y*.

.. function:: void acb_mul_onei(acb_t z, const acb_t x)

    Sets *z* to *x* multiplied by the imaginary unit.

.. function:: void acb_div_onei(acb_t z, const acb_t x)

    Sets *z* to *x* divided by the imaginary unit.

.. function:: void acb_mul_ui(acb_t z, const acb_t x, ulong y, slong prec)

.. function:: void acb_mul_si(acb_t z, const acb_t x, slong y, slong prec)

.. function:: void acb_mul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)

.. function:: void acb_mul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

    Sets *z* to the product of *x* and *y*.

.. function:: void acb_mul(acb_t z, const acb_t x, const acb_t y, slong prec)

    Sets *z* to the product of *x* and *y*. If at least one part of
    *x* or *y* is zero, the operations is reduced to two real multiplications.
    If *x* and *y* are the same pointers, they are assumed to represent
    the same mathematical quantity and the squaring formula is used.

.. function:: void acb_mul_2exp_si(acb_t z, const acb_t x, slong e)

.. function:: void acb_mul_2exp_fmpz(acb_t z, const acb_t x, const fmpz_t e)

    Sets *z* to *x* multiplied by `2^e`, without rounding.

.. function:: void acb_sqr(acb_t z, const acb_t x, slong prec)

    Sets *z* to *x* squared.

.. function:: void acb_cube(acb_t z, const acb_t x, slong prec)

    Sets *z* to *x* cubed, computed efficiently using two real squarings,
    two real multiplications, and scalar operations.

.. function:: void acb_addmul(acb_t z, const acb_t x, const acb_t y, slong prec)

.. function:: void acb_addmul_ui(acb_t z, const acb_t x, ulong y, slong prec)

.. function:: void acb_addmul_si(acb_t z, const acb_t x, slong y, slong prec)

.. function:: void acb_addmul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)

.. function:: void acb_addmul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

    Sets *z* to *z* plus the product of *x* and *y*.

.. function:: void acb_submul(acb_t z, const acb_t x, const acb_t y, slong prec)

.. function:: void acb_submul_ui(acb_t z, const acb_t x, ulong y, slong prec)

.. function:: void acb_submul_si(acb_t z, const acb_t x, slong y, slong prec)

.. function:: void acb_submul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)

.. function:: void acb_submul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

    Sets *z* to *z* minus the product of *x* and *y*.

.. function:: void acb_inv(acb_t z, const acb_t x, slong prec)

    Sets *z* to the multiplicative inverse of *x*.

.. function:: void acb_div_ui(acb_t z, const acb_t x, ulong y, slong prec)

.. function:: void acb_div_si(acb_t z, const acb_t x, slong y, slong prec)

.. function:: void acb_div_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)

.. function:: void acb_div_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

.. function:: void acb_div(acb_t z, const acb_t x, const acb_t y, slong prec)

    Sets *z* to the quotient of *x* and *y*.

Dot product
-------------------------------------------------------------------------------

.. function:: void acb_dot_precise(acb_t res, const acb_t s, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)
              void acb_dot_simple(acb_t res, const acb_t s, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)
              void acb_dot(acb_t res, const acb_t s, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)

    Computes the dot product of the vectors *x* and *y*, setting
    *res* to `s + (-1)^{subtract} \sum_{i=0}^{len-1} x_i y_i`.

    The initial term *s* is optional and can be
    omitted by passing *NULL* (equivalently, `s = 0`).
    The parameter *subtract* must be 0 or 1.
    The length *len* is allowed to be negative, which is equivalent
    to a length of zero.
    The parameters *xstep* or *ystep* specify a step length for
    traversing subsequences of the vectors *x* and *y*; either can be
    negative to step in the reverse direction starting from
    the initial pointer.
    Aliasing is allowed between *res* and *s* but not between
    *res* and the entries of *x* and *y*.

    The default version determines the optimal precision for each term
    and performs all internal calculations using mpn arithmetic
    with minimal overhead. This is the preferred way to compute a
    dot product; it is generally much faster and more precise
    than a simple loop.

    The *simple* version performs fused multiply-add operations in
    a simple loop. This can be used for
    testing purposes and is also used as a fallback by the
    default version when the exponents are out of range
    for the optimized code.

    The *precise* version computes the dot product exactly up to the
    final rounding. This can be extremely slow and is only intended
    for testing.

.. function:: void acb_approx_dot(acb_t res, const acb_t s, int subtract, acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)

    Computes an approximate dot product *without error bounds*.
    The radii of the inputs are ignored (only the midpoints are read)
    and only the midpoint of the output is written.

.. function:: void acb_dot_ui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
              void acb_dot_si(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const slong * y, slong ystep, slong len, slong prec)
              void acb_dot_uiui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
              void acb_dot_siui(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
              void acb_dot_fmpz(acb_t res, const acb_t initial, int subtract, acb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec)

    Equivalent to :func:`acb_dot`, but with integers in the array *y*.
    The *uiui* and *siui* versions take an array of double-limb integers
    as input; the *siui* version assumes that these represent signed
    integers in two's complement form.

Mathematical constants
-------------------------------------------------------------------------------

.. function:: void acb_const_pi(acb_t y, slong prec)

    Sets *y* to the constant `\pi`.

Powers and roots
-------------------------------------------------------------------------------

.. function:: void acb_sqrt(acb_t r, const acb_t z, slong prec)

    Sets *r* to the square root of *z*.
    If either the real or imaginary part is exactly zero, only
    a single real square root is needed. Generally, we use the formula
    `\sqrt{a+bi} = u/2 + ib/u, u = \sqrt{2(|a+bi|+a)}`,
    requiring two real square root extractions.

.. function:: void acb_sqrt_analytic(acb_t r, const acb_t z, int analytic, slong prec)

    Computes the square root. If *analytic* is set, gives a NaN-containing
    result if *z* touches the branch cut.

.. function:: void acb_rsqrt(acb_t r, const acb_t z, slong prec)

    Sets *r* to the reciprocal square root of *z*.
    If either the real or imaginary part is exactly zero, only
    a single real reciprocal square root is needed. Generally, we use the
    formula `1/\sqrt{a+bi} = ((a+r) - bi)/v, r = |a+bi|, v = \sqrt{r |a+bi+r|^2}`,
    requiring one real square root and one real reciprocal square root.

.. function:: void acb_rsqrt_analytic(acb_t r, const acb_t z, int analytic, slong prec)

    Computes the reciprocal square root. If *analytic* is set, gives a
    NaN-containing result if *z* touches the branch cut.

.. function:: void acb_quadratic_roots_fmpz(acb_t r1, acb_t r2, const fmpz_t a, const fmpz_t b, const fmpz_t c, slong prec)

    Sets *r1* and *r2* to the roots of the quadratic polynomial
    `ax^2 + bx + c`. Requires that *a* is nonzero.
    This function is implemented so that both roots are computed accurately
    even when direct use of the quadratic formula would lose accuracy.

.. function:: void acb_root_ui(acb_t r, const acb_t z, ulong k, slong prec)

    Sets *r* to the principal *k*-th root of *z*.

.. function:: void acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, slong prec)

.. function:: void acb_pow_ui(acb_t y, const acb_t b, ulong e, slong prec)

.. function:: void acb_pow_si(acb_t y, const acb_t b, slong e, slong prec)

    Sets `y = b^e` using binary exponentiation (with an initial division
    if `e < 0`). Note that these functions can get slow if the exponent is
    extremely large (in such cases :func:`acb_pow` may be superior).

.. function:: void acb_pow_arb(acb_t z, const acb_t x, const arb_t y, slong prec)

.. function:: void acb_pow(acb_t z, const acb_t x, const acb_t y, slong prec)

    Sets `z = x^y`, computed using binary exponentiation if `y` if
    a small exact integer, as `z = (x^{1/2})^{2y}` if `y` is a small exact
    half-integer, and generally as `z = \exp(y \log x)`.

.. function:: void acb_pow_analytic(acb_t r, const acb_t x, const acb_t y, int analytic, slong prec)

    Computes the power `x^y`. If *analytic* is set, gives a
    NaN-containing result if *x* touches the branch cut (unless *y* is
    an integer).

.. function:: void acb_unit_root(acb_t res, ulong order, slong prec)

    Sets *res* to `\exp(\frac{2i\pi}{\mathrm{order}})` to precision *prec*.

Exponentials and logarithms
-------------------------------------------------------------------------------

.. function:: void acb_exp(acb_t y, const acb_t z, slong prec)

    Sets *y* to the exponential function of *z*, computed as
    `\exp(a+bi) = \exp(a) \left( \cos(b) + \sin(b) i \right)`.

.. function:: void acb_exp_pi_i(acb_t y, const acb_t z, slong prec)

    Sets *y* to `\exp(\pi i z)`.

.. function:: void acb_exp_invexp(acb_t s, acb_t t, const acb_t z, slong prec)

    Sets `s = \exp(z)` and `t = \exp(-z)`.

.. function:: void acb_expm1(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\exp(z)-1`, using a more accurate method when `z \approx 0`.

.. function:: void acb_log(acb_t y, const acb_t z, slong prec)

    Sets *y* to the principal branch of the natural logarithm of *z*,
    computed as
    `\log(a+bi) = \frac{1}{2} \log(a^2 + b^2) + i \operatorname{arg}(a+bi)`.

.. function:: void acb_log_analytic(acb_t r, const acb_t z, int analytic, slong prec)

    Computes the natural logarithm. If *analytic* is set, gives a
    NaN-containing result if *z* touches the branch cut.

.. function:: void acb_log1p(acb_t z, const acb_t x, slong prec)

    Sets `z = \log(1+x)`, computed accurately when `x \approx 0`.

Trigonometric functions
-------------------------------------------------------------------------------

.. function:: void acb_sin(acb_t s, const acb_t z, slong prec)

.. function:: void acb_cos(acb_t c, const acb_t z, slong prec)

.. function:: void acb_sin_cos(acb_t s, acb_t c, const acb_t z, slong prec)

    Sets `s = \sin(z)`, `c = \cos(z)`, evaluated as
    `\sin(a+bi) = \sin(a)\cosh(b) + i \cos(a)\sinh(b)`,
    `\cos(a+bi) = \cos(a)\cosh(b) - i \sin(a)\sinh(b)`.

.. function:: void acb_tan(acb_t s, const acb_t z, slong prec)

    Sets `s = \tan(z) = \sin(z) / \cos(z)`. For large imaginary parts,
    the function is evaluated in a numerically stable way as `\pm i`
    plus a decreasing exponential factor.

.. function:: void acb_cot(acb_t s, const acb_t z, slong prec)

    Sets `s = \cot(z) = \cos(z) / \sin(z)`. For large imaginary parts,
    the function is evaluated in a numerically stable way as `\pm i`
    plus a decreasing exponential factor.

.. function:: void acb_sin_pi(acb_t s, const acb_t z, slong prec)

.. function:: void acb_cos_pi(acb_t s, const acb_t z, slong prec)

.. function:: void acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, slong prec)

    Sets `s = \sin(\pi z)`, `c = \cos(\pi z)`, evaluating the trigonometric
    factors of the real and imaginary part accurately via :func:`arb_sin_cos_pi`.

.. function:: void acb_tan_pi(acb_t s, const acb_t z, slong prec)

    Sets `s = \tan(\pi z)`. Uses the same algorithm as :func:`acb_tan`,
    but evaluates the sine and cosine accurately via :func:`arb_sin_cos_pi`.

.. function:: void acb_cot_pi(acb_t s, const acb_t z, slong prec)

    Sets `s = \cot(\pi z)`. Uses the same algorithm as :func:`acb_cot`,
    but evaluates the sine and cosine accurately via :func:`arb_sin_cos_pi`.

.. function:: void acb_sec(acb_t res, const acb_t z, slong prec)

    Computes `\sec(z) = 1 / \cos(z)`.

.. function:: void acb_csc(acb_t res, const acb_t z, slong prec)

    Computes `\csc(x) = 1 / \sin(z)`.

.. function:: void acb_csc_pi(acb_t res, const acb_t z, slong prec)

    Computes `\csc(\pi x) = 1 / \sin(\pi z)`. Evaluates the sine accurately
    via :func:`acb_sin_pi`.

.. function:: void acb_sinc(acb_t s, const acb_t z, slong prec)

    Sets `s = \operatorname{sinc}(x) = \sin(z) / z`.

.. function:: void acb_sinc_pi(acb_t s, const acb_t z, slong prec)

    Sets `s = \operatorname{sinc}(\pi x) = \sin(\pi z) / (\pi z)`.

Inverse trigonometric functions
-------------------------------------------------------------------------------

.. function:: void acb_asin(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\operatorname{asin}(z) = -i \log(iz + \sqrt{1-z^2})`.

.. function:: void acb_acos(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\operatorname{acos}(z) = \tfrac{1}{2} \pi - \operatorname{asin}(z)`.

.. function:: void acb_atan(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\operatorname{atan}(z) = \tfrac{1}{2} i (\log(1-iz)-\log(1+iz))`.

Hyperbolic functions
-------------------------------------------------------------------------------

.. function:: void acb_sinh(acb_t s, const acb_t z, slong prec)

.. function:: void acb_cosh(acb_t c, const acb_t z, slong prec)

.. function:: void acb_sinh_cosh(acb_t s, acb_t c, const acb_t z, slong prec)

.. function:: void acb_tanh(acb_t s, const acb_t z, slong prec)

.. function:: void acb_coth(acb_t s, const acb_t z, slong prec)

    Respectively computes `\sinh(z) = -i\sin(iz)`, `\cosh(z) = \cos(iz)`,
    `\tanh(z) = -i\tan(iz)`, `\coth(z) = i\cot(iz)`.

.. function:: void acb_sech(acb_t res, const acb_t z, slong prec)

    Computes `\operatorname{sech}(z) = 1 / \cosh(z)`.

.. function:: void acb_csch(acb_t res, const acb_t z, slong prec)

    Computes `\operatorname{csch}(z) = 1 / \sinh(z)`.

Inverse hyperbolic functions
-------------------------------------------------------------------------------

.. function:: void acb_asinh(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\operatorname{asinh}(z) = -i \operatorname{asin}(iz)`.

.. function:: void acb_acosh(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\operatorname{acosh}(z) = \log(z + \sqrt{z+1} \sqrt{z-1})`.

.. function:: void acb_atanh(acb_t res, const acb_t z, slong prec)

    Sets *res* to `\operatorname{atanh}(z) = -i \operatorname{atan}(iz)`.

Lambert W function
-------------------------------------------------------------------------------

.. function:: void acb_lambertw_asymp(acb_t res, const acb_t z, const fmpz_t k, slong L, slong M, slong prec)

    Sets *res* to the Lambert W function `W_k(z)` computed using *L* and *M*
    terms in the bivariate series giving the asymptotic expansion at
    zero or infinity. This algorithm is valid
    everywhere, but the error bound is only finite when `|\log(z)|` is
    sufficiently large.

.. function:: int acb_lambertw_check_branch(const acb_t w, const fmpz_t k, slong prec)

    Tests if *w* definitely lies in the image of the branch `W_k(z)`.
    This function is used internally to verify that a computed approximation
    of the Lambert W function lies on the intended branch. Note that this will
    necessarily evaluate to false for points exactly on (or overlapping) the
    branch cuts, where a different algorithm has to be used.

.. function:: void acb_lambertw_bound_deriv(mag_t res, const acb_t z, const acb_t ez1, const fmpz_t k)

    Sets *res* to an upper bound for `|W_k'(z)|`. The input *ez1* should
    contain the precomputed value of `ez+1`.

    Along the real line, the directional derivative of `W_k(z)` is understood
    to be taken. As a result, the user must handle the branch cut
    discontinuity separately when using this function to bound perturbations
    in the value of `W_k(z)`.

.. function:: void acb_lambertw(acb_t res, const acb_t z, const fmpz_t k, int flags, slong prec)

    Sets *res* to the Lambert W function `W_k(z)` where the index *k* selects
    the branch (with `k = 0` giving the principal branch).
    The placement of branch cuts follows [CGHJK1996]_.

    If *flags* is nonzero, nonstandard branch cuts are used.

    If *flags* is set to *ACB_LAMBERTW_LEFT*, computes `W_{\mathrm{left}|k}(z)`
    which corresponds to `W_k(z)` in the upper
    half plane and `W_{k+1}(z)` in the lower half plane, connected continuously
    to the left of the branch points.
    In other words, the branch cut on `(-\infty,0)` is rotated counterclockwise
    to `(0,+\infty)`.
    (For `k = -1` and `k = 0`, there is also a branch cut on `(-1/e,0)`,
    continuous from below instead of from above to maintain counterclockwise
    continuity.)

    If *flags* is set to *ACB_LAMBERTW_MIDDLE*, computes
    `W_{\mathrm{middle}}(z)` which corresponds to
    `W_{-1}(z)` in the upper half plane and `W_{1}(z)` in the lower half
    plane, connected continuously through `(-1/e,0)` with branch cuts
    on `(-\infty,-1/e)` and `(0,+\infty)`. `W_{\mathrm{middle}}(z)` extends the
    real analytic function `W_{-1}(x)` defined on `(-1/e,0)` to a complex
    analytic function, whereas the standard branch `W_{-1}(z)` has a branch
    cut along the real segment.

    The algorithm used to compute the Lambert W function is described
    in [Joh2017b]_.

Rising factorials
-------------------------------------------------------------------------------

.. function:: void acb_rising_ui(acb_t z, const acb_t x, ulong n, slong prec)
              void acb_rising(acb_t z, const acb_t x, const acb_t n, slong prec)

    Computes the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)`.
    These functions are aliases for :func:`acb_hypgeom_rising_ui`
    and :func:`acb_hypgeom_rising`.

.. function :: void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)

    Letting `u(x) = x (x+1) (x+2) \cdots (x+n-1)`, simultaneously compute
    `u(x)` and `v(x) = u'(x)`.
    This function is a wrapper of :func:`acb_hypgeom_rising_ui_jet`.

.. function :: void acb_rising_ui_get_mag(mag_t bound, const acb_t x, ulong n)

    Computes an upper bound for the absolute value of
    the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)`.
    Not currently optimized for large *n*.

Gamma function
-------------------------------------------------------------------------------

.. function:: void acb_gamma(acb_t y, const acb_t x, slong prec)

    Computes the gamma function `y = \Gamma(x)`.
    This is an alias for :func:`acb_hypgeom_gamma`.

.. function:: void acb_rgamma(acb_t y, const acb_t x, slong prec)

    Computes the reciprocal gamma function  `y = 1/\Gamma(x)`,
    avoiding division by zero at the poles of the gamma function.
    This is an alias for :func:`acb_hypgeom_rgamma`.

.. function:: void acb_lgamma(acb_t y, const acb_t x, slong prec)

    Computes the logarithmic gamma function `y = \log \Gamma(x)`.
    This is an alias for :func:`acb_hypgeom_lgamma`.

    The branch cut of the logarithmic gamma function is placed on the
    negative half-axis, which means that
    `\log \Gamma(z) + \log z = \log \Gamma(z+1)` holds for all `z`,
    whereas `\log \Gamma(z) \ne \log(\Gamma(z))` in general.
    In the left half plane, the reflection formula with correct
    branch structure is evaluated via :func:`acb_log_sin_pi`.

.. function:: void acb_digamma(acb_t y, const acb_t x, slong prec)

    Computes the digamma function `y = \psi(x) = (\log \Gamma(x))' = \Gamma'(x) / \Gamma(x)`.

.. function:: void acb_log_sin_pi(acb_t res, const acb_t z, slong prec)

    Computes the logarithmic sine function defined by

    .. math ::

        S(z) = \log(\pi) - \log \Gamma(z) + \log \Gamma(1-z)

    which is equal to

    .. math ::

        S(z) = \int_{1/2}^z \pi \cot(\pi t) dt

    where the path of integration goes through the upper half plane
    if `0 < \arg(z) \le \pi` and through the lower half plane
    if `-\pi < \arg(z) \le 0`. Equivalently,

    .. math ::

        S(z) = \log(\sin(\pi(z-n))) \mp n \pi i, \quad n = \lfloor \operatorname{re}(z) \rfloor

    where the negative sign is taken if `0 < \arg(z) \le \pi`
    and the positive sign is taken otherwise (if the interval `\arg(z)`
    does not certainly satisfy either condition, the union of
    both cases is computed).
    After subtracting *n*, we have `0 \le \operatorname{re}(z) < 1`. In
    this strip, we use
    use `S(z) = \log(\sin(\pi(z)))` if the imaginary part of *z* is small.
    Otherwise, we use `S(z) = i \pi (z-1/2) + \log((1+e^{-2i\pi z})/2)`
    in the lower half-plane and the conjugated expression in the upper
    half-plane to avoid exponent overflow.

    The function is evaluated at the midpoint and the propagated error
    is computed from `S'(z)` to get a continuous change
    when `z` is non-real and `n` spans more than one possible integer value.

.. function:: void acb_polygamma(acb_t res, const acb_t s, const acb_t z, slong prec)

    Sets *res* to the value of the generalized polygamma function `\psi(s,z)`.

    If *s* is a nonnegative order, this is simply the *s*-order derivative
    of the digamma function. If `s = 0`, this function simply
    calls the digamma function internally. For integers `s \ge 1`,
    it calls the Hurwitz zeta function. Note that for small integers
    `s \ge 1`, it can be faster to use
    :func:`acb_poly_digamma_series` and read off the coefficients.

    The generalization to other values of *s* is due to
    Espinosa and Moll [EM2004]_:

    .. math ::

        \psi(s,z) = \frac{\zeta'(s+1,z) + (\gamma + \psi(-s)) \zeta(s+1,z)}{\Gamma(-s)}

.. function:: void acb_barnes_g(acb_t res, const acb_t z, slong prec)

.. function:: void acb_log_barnes_g(acb_t res, const acb_t z, slong prec)

    Computes Barnes *G*-function or the logarithmic Barnes *G*-function,
    respectively. The logarithmic version has branch cuts on the negative
    real axis and is continuous elsewhere in the complex plane,
    in analogy with the logarithmic gamma function. The functional
    equation

    .. math ::

        \log G(z+1) = \log \Gamma(z) + \log G(z).

    holds for all *z*.

    For small integers, we directly use the recurrence
    relation `G(z+1) = \Gamma(z) G(z)` together with the initial value
    `G(1) = 1`. For general *z*, we use the formula

    .. math ::

        \log G(z) = (z-1) \log \Gamma(z) - \zeta'(-1,z) + \zeta'(-1).

Zeta function
-------------------------------------------------------------------------------

.. function:: void acb_zeta(acb_t z, const acb_t s, slong prec)

    Sets *z* to the value of the Riemann zeta function `\zeta(s)`.
    Note: for computing derivatives with respect to `s`,
    use :func:`acb_poly_zeta_series` or related methods.

    This is a wrapper of :func:`acb_dirichlet_zeta`.

.. function:: void acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec)

    Sets *z* to the value of the Hurwitz zeta function `\zeta(s, a)`.
    Note: for computing derivatives with respect to `s`,
    use :func:`acb_poly_zeta_series` or related methods.

    This is a wrapper of :func:`acb_dirichlet_hurwitz`.

.. function:: void acb_bernoulli(acb_t z, const acb_t s, slong prec)

    Sets *z* to the value of the Bernoulli function `B(s)`.

    This is a wrapper of :func:`acb_dirichlet_bernoulli`.

.. function:: void acb_bernoulli_gen(acb_t z, const acb_t s, const acb_t a, slong prec)

    Sets *z* to the value of the generalized Bernoulli function `B(s, a)`.

    This is a wrapper of :func:`acb_dirichlet_bernoulli_gen`.

.. function:: void acb_bernoulli_poly_ui(acb_t res, ulong n, const acb_t x, slong prec)

    Sets *res* to the value of the Bernoulli polynomial `B_n(x)`.

    Warning: this function is only fast if either *n* or *x* is a small integer.

    This function reads Bernoulli numbers from the global cache if they
    are already cached, but does not automatically extend the cache by itself.

Polylogarithms
-------------------------------------------------------------------------------

.. function:: void acb_polylog(acb_t w, const acb_t s, const acb_t z, slong prec)

.. function:: void acb_polylog_si(acb_t w, slong s, const acb_t z, slong prec)

    Sets *w* to the polylogarithm `\operatorname{Li}_s(z)`.

Arithmetic-geometric mean
-------------------------------------------------------------------------------

See :ref:`algorithms_agm` for implementation details.

.. function:: void acb_agm1(acb_t m, const acb_t z, slong prec)

    Sets *m* to the arithmetic-geometric mean `M(z) = \operatorname{agm}(1,z)`,
    defined such that the function is continuous in the complex plane except for
    a branch cut along the negative half axis (where it is continuous
    from above). This corresponds to always choosing an "optimal" branch for
    the square root in the arithmetic-geometric mean iteration.

.. function:: void acb_agm1_cpx(acb_ptr m, const acb_t z, slong len, slong prec)

    Sets the coefficients in the array *m* to the power series expansion of the
    arithmetic-geometric mean at the point *z* truncated to length *len*, i.e.
    `M(z+x) \in \mathbb{C}[[x]]`.

.. function:: void acb_agm(acb_t m, const acb_t x, const acb_t y, slong prec)

    Sets *m* to the arithmetic-geometric mean of *x* and *y*. The square
    roots in the AGM iteration are chosen so as to form the "optimal"
    AGM sequence. This gives a well-defined function of *x* and *y* except
    when `x / y` is a negative real number, in which case there are two
    optimal AGM sequences. In that case, an arbitrary but consistent
    choice is made (if a decision cannot be made due to inexact arithmetic,
    the union of both choices is returned).

Other special functions
-------------------------------------------------------------------------------

.. function:: void acb_chebyshev_t_ui(acb_t a, ulong n, const acb_t x, slong prec)

.. function:: void acb_chebyshev_u_ui(acb_t a, ulong n, const acb_t x, slong prec)

    Evaluates the Chebyshev polynomial of the first kind `a = T_n(x)`
    or the Chebyshev polynomial of the second kind `a = U_n(x)`.

.. function:: void acb_chebyshev_t2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec)

.. function:: void acb_chebyshev_u2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec)

    Simultaneously evaluates `a = T_n(x), b = T_{n-1}(x)` or
    `a = U_n(x), b = U_{n-1}(x)`.
    Aliasing between *a*, *b* and *x* is not permitted.

Piecewise real functions
-------------------------------------------------------------------------------

The following methods extend common piecewise real functions to piecewise
complex analytic functions, useful together with the
:ref:`acb_calc.h <acb-calc>` module.
If *analytic* is set, evaluation on a discontinuity or non-analytic point
gives a NaN result.

.. function:: void acb_real_abs(acb_t res, const acb_t z, int analytic, slong prec)

    The absolute value is extended to `+z` in the right half plane and
    `-z` in the left half plane, with a discontinuity on the vertical line
    `\operatorname{Re}(z) = 0`.

.. function:: void acb_real_sgn(acb_t res, const acb_t z, int analytic, slong prec)

    The sign function is extended to `+1` in the right half plane and
    `-1` in the left half plane, with a discontinuity on the vertical line
    `\operatorname{Re}(z) = 0`.
    If *analytic* is not set, this is effectively the same function as
    :func:`acb_csgn`.

.. function:: void acb_real_heaviside(acb_t res, const acb_t z, int analytic, slong prec)

    The Heaviside step function (or unit step function) is extended to `+1` in
    the right half plane and `0` in the left half plane, with a discontinuity on
    the vertical line `\operatorname{Re}(z) = 0`.

.. function:: void acb_real_floor(acb_t res, const acb_t z, int analytic, slong prec)

    The floor function is extended to a piecewise constant function
    equal to `n` in the strips with real part `(n,n+1)`, with discontinuities
    on the vertical lines `\operatorname{Re}(z) = n`.

.. function:: void acb_real_ceil(acb_t res, const acb_t z, int analytic, slong prec)

    The ceiling function is extended to a piecewise constant function
    equal to `n+1` in the strips with real part `(n,n+1)`, with discontinuities
    on the vertical lines `\operatorname{Re}(z) = n`.

.. function:: void acb_real_max(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec)

    The real function `\max(x,y)` is extended to a piecewise analytic function
    of two variables by returning `x` when
    `\operatorname{Re}(x) \ge \operatorname{Re}(y)`
    and returning `y` when `\operatorname{Re}(x) < \operatorname{Re}(y)`,
    with discontinuities where `\operatorname{Re}(x) = \operatorname{Re}(y)`.

.. function:: void acb_real_min(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec)

    The real function `\min(x,y)` is extended to a piecewise analytic function
    of two variables by returning `x` when
    `\operatorname{Re}(x) \le \operatorname{Re}(y)`
    and returning `y` when `\operatorname{Re}(x) > \operatorname{Re}(y)`,
    with discontinuities where `\operatorname{Re}(x) = \operatorname{Re}(y)`.

.. function:: void acb_real_sqrtpos(acb_t res, const acb_t z, int analytic, slong prec)

    Extends the real square root function on `[0,+\infty)` to the usual
    complex square root on the cut plane. Like :func:`arb_sqrtpos`, only
    the nonnegative part of *z* is considered if *z* is purely real
    and *analytic* is not set. This is useful for integrating `\sqrt{f(x)}`
    where it is known that `f(x) \ge 0`: unlike :func:`acb_sqrt_analytic`,
    no spurious imaginary terms `[\pm \varepsilon] i` are created when the
    balls computed for `f(x)` straddle zero.

Vector functions
-------------------------------------------------------------------------------

.. function:: void _acb_vec_zero(acb_ptr A, slong n)

    Sets all entries in *vec* to zero.

.. function:: int _acb_vec_is_zero(acb_srcptr vec, slong len)

    Returns nonzero iff all entries in *x* are zero.

.. function:: int _acb_vec_is_real(acb_srcptr v, slong len)

    Returns nonzero iff all entries in *x* have zero imaginary part.

.. function:: void _acb_vec_set(acb_ptr res, acb_srcptr vec, slong len)

    Sets *res* to a copy of *vec*.

.. function:: void _acb_vec_set_round(acb_ptr res, acb_srcptr vec, slong len, slong prec)

    Sets *res* to a copy of *vec*, rounding each entry to *prec* bits.

.. function:: void _acb_vec_swap(acb_ptr vec1, acb_ptr vec2, slong len)

    Swaps the entries of *vec1* and *vec2*.

.. function:: void _acb_vec_neg(acb_ptr res, acb_srcptr vec, slong len)

.. function:: void _acb_vec_add(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, slong len, slong prec)

.. function:: void _acb_vec_sub(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, slong len, slong prec)

.. function:: void _acb_vec_scalar_submul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)

.. function:: void _acb_vec_scalar_addmul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)

.. function:: void _acb_vec_scalar_mul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)

.. function:: void _acb_vec_scalar_mul_ui(acb_ptr res, acb_srcptr vec, slong len, ulong c, slong prec)

.. function:: void _acb_vec_scalar_mul_2exp_si(acb_ptr res, acb_srcptr vec, slong len, slong c)

.. function:: void _acb_vec_scalar_mul_onei(acb_ptr res, acb_srcptr vec, slong len)

.. function:: void _acb_vec_scalar_div_ui(acb_ptr res, acb_srcptr vec, slong len, ulong c, slong prec)

.. function:: void _acb_vec_scalar_div(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)

.. function:: void _acb_vec_scalar_mul_arb(acb_ptr res, acb_srcptr vec, slong len, const arb_t c, slong prec)

.. function:: void _acb_vec_scalar_div_arb(acb_ptr res, acb_srcptr vec, slong len, const arb_t c, slong prec)

.. function:: void _acb_vec_scalar_mul_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec)

.. function:: void _acb_vec_scalar_div_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec)

   Performs the respective scalar operation elementwise.

.. function:: slong _acb_vec_bits(acb_srcptr vec, slong len)

    Returns the maximum of :func:`arb_bits` for all entries in *vec*.

.. function:: void _acb_vec_set_powers(acb_ptr xs, const acb_t x, slong len, slong prec)

    Sets *xs* to the powers `1, x, x^2, \ldots, x^{len-1}`.

.. function:: void _acb_vec_unit_roots(acb_ptr z, slong order, slong len, slong prec)

    Sets *z* to the powers `1,z,z^2,\dots z^{\mathrm{len}-1}` where `z=\exp(\frac{2i\pi}{\mathrm{order}})` to precision *prec*.
    *order* can be taken negative.

    In order to avoid precision loss, this function does not simply compute powers of a primitive root.

.. function:: void _acb_vec_add_error_arf_vec(acb_ptr res, arf_srcptr err, slong len)

.. function:: void _acb_vec_add_error_mag_vec(acb_ptr res, mag_srcptr err, slong len)

    Adds the magnitude of each entry in *err* to the radius of the
    corresponding entry in *res*.

.. function:: void _acb_vec_indeterminate(acb_ptr vec, slong len)

    Applies :func:`acb_indeterminate` elementwise.

.. function:: void _acb_vec_trim(acb_ptr res, acb_srcptr vec, slong len)

    Applies :func:`acb_trim` elementwise.

.. function:: int _acb_vec_get_unique_fmpz_vec(fmpz * res,  acb_srcptr vec, slong len)

    Calls :func:`acb_get_unique_fmpz` elementwise and returns nonzero if
    all entries can be rounded uniquely to integers. If any entry in *vec*
    cannot be rounded uniquely to an integer, returns zero.

.. function:: void _acb_vec_sort_pretty(acb_ptr vec, slong len)

    Sorts the vector of complex numbers based on the real and imaginary parts.
    This is intended to reveal structure when printing a set of complex numbers,
    not to apply an order relation in a rigorous way.

