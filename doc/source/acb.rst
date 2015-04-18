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

.. function:: void acb_init(arb_t x)

    Initializes the variable *x* for use, and sets its value to zero.

.. function:: void acb_clear(acb_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: acb_ptr _acb_vec_init(long n)

    Returns a pointer to an array of *n* initialized *acb_struct*:s.

.. function:: void _acb_vec_clear(acb_ptr v, long n)

    Clears an array of *n* initialized *acb_struct*:s.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: int acb_is_zero(const acb_t z)

    Returns nonzero iff *z* is zero.

.. function:: int acb_is_one(const acb_t z)

    Returns nonzero iff *z* is exactly 1.

.. function:: int acb_is_exact(const acb_t z)

    Returns nonzero iff *z* is exact.

.. function:: int acb_is_int(const acb_t z)

    Returns nonzero iff *z* is an exact integer.

.. function:: void acb_zero(acb_t z)

.. function:: void acb_one(acb_t z)

.. function:: void acb_onei(acb_t z)

    Sets *z* respectively to 0, 1, `i = \sqrt{-1}`.

.. function:: void acb_set(acb_t z, const acb_t x)

.. function:: void acb_set_ui(acb_t z, long x)

.. function:: void acb_set_si(acb_t z, long x)

.. function:: void acb_set_fmpz(acb_t z, const fmpz_t x)

.. function:: void acb_set_arb(acb_t z, const arb_t c)

    Sets *z* to the value of *x*.

.. function:: void acb_set_fmpq(acb_t z, const fmpq_t x, long prec)

.. function:: void acb_set_round(acb_t z, const acb_t x, long prec)

.. function:: void acb_set_round_fmpz(acb_t z, const fmpz_t x, long prec)

.. function:: void acb_set_round_arb(acb_t z, const arb_t x, long prec)

    Sets *z* to *x*, rounded to *prec* bits.

.. function:: void acb_swap(acb_t z, acb_t x)

    Swaps *z* and *x* efficiently.


Input and output
-------------------------------------------------------------------------------

.. function:: void acb_print(const acb_t x)

    Prints the internal representation of *x*.

.. function:: void acb_printd(const acb_t z, long digits)

    Prints *x* in decimal. The printed value of the radius is not adjusted
    to compensate for the fact that the binary-to-decimal conversion
    of both the midpoint and the radius introduces additional error.


Random number generation
-------------------------------------------------------------------------------

.. function:: void acb_randtest(acb_t z, flint_rand_t state, long prec, long mag_bits)

    Generates a random complex number by generating separate random
    real and imaginary parts.

.. function:: void acb_randtest_special(acb_t z, flint_rand_t state, long prec, long mag_bits)

    Generates a random complex number by generating separate random
    real and imaginary parts. Also generates NaNs and infinities.

.. function:: void acb_randtest_precise(acb_t z, flint_rand_t state, long prec, long mag_bits)

    Generates a random complex number with precise real and imaginary parts.

Precision and comparisons
-------------------------------------------------------------------------------

.. function:: int acb_equal(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* are identical as sets, i.e.
    if the real and imaginary parts are equal as balls.

    Note that this is not the same thing as testing whether both
    *x* and *y* certainly represent the same complex number, unless
    either *x* or *y* is exact (and neither contains NaN).
    To test whether both operands *might* represent the same mathematical
    quantity, use :func:`acb_overlaps` or :func:`acb_contains`,
    depending on the circumstance.

.. function:: int acb_eq(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* are certainly equal, as determined
    by testing that :func:`arb_eq` holds for both the real and imaginary
    parts.

.. function:: int acb_ne(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* are certainly not equal, as determined
    by testing that :func:`arb_ne` holds for either the real or imaginary parts.

.. function:: int acb_overlaps(const acb_t x, const acb_t y)

    Returns nonzero iff *x* and *y* have some point in common.

.. function:: void acb_get_abs_ubound_arf(arf_t u, const acb_t z, long prec)

    Sets *u* to an upper bound for the absolute value of *z*, computed
    using a working precision of *prec* bits.

.. function:: void acb_get_abs_lbound_arf(arf_t u, const acb_t z, long prec)

    Sets *u* to a lower bound for the absolute value of *z*, computed
    using a working precision of *prec* bits.

.. function:: void acb_get_rad_ubound_arf(arf_t u, const acb_t z, long prec)

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

.. function:: long acb_rel_error_bits(const acb_t x)

    Returns the effective relative error of *x* measured in bits.
    This is computed as if calling :func:`arb_rel_error_bits` on the
    real ball whose midpoint is the larger out of the real and imaginary
    midpoints of *x*, and whose radius is the larger out of the real
    and imaginary radiuses of *x*.

.. function:: long acb_rel_accuracy_bits(const arb_t x)

    Returns the effective relative accuracy of *x* measured in bits,
    equal to the negative of the return value from :func:`acb_rel_error_bits`.

.. function:: long acb_bits(const acb_t x)

    Returns the maximum of *arb_bits* applied to the real
    and imaginary parts of *x*, i.e. the minimum precision sufficient
    to represent *x* exactly.

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

.. function:: void acb_arg(arb_t r, const acb_t z, long prec)

    Sets *r* to a real interval containing the complex argument (phase) of *z*.
    We define the complex argument have a discontinuity on `(-\infty,0]`, with
    the special value `\operatorname{arg}(0) = 0`, and
    `\operatorname{arg}(a+0i) = \pi` for `a < 0`. Equivalently, if
    `z = a+bi`, the argument is given by `\operatorname{atan2}(b,a)`
    (see :func:`arb_atan2`).

.. function:: void acb_abs(arb_t r, const acb_t z, long prec)

    Sets *r* to the absolute value of *z*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void acb_neg(acb_t z, const acb_t x)

    Sets *z* to the negation of *x*.

.. function:: void acb_conj(acb_t z, const acb_t x)

    Sets *z* to the complex conjugate of *x*.

.. function:: void acb_add_ui(acb_t z, const acb_t x, ulong y, long prec)

.. function:: void acb_add_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)

.. function:: void acb_add_arb(acb_t z, const acb_t x, const arb_t y, long prec)

.. function:: void acb_add(acb_t z, const acb_t x, const acb_t y, long prec)

    Sets *z* to the sum of *x* and *y*.

.. function:: void acb_sub_ui(acb_t z, const acb_t x, ulong y, long prec)

.. function:: void acb_sub_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)

.. function:: void acb_sub_arb(acb_t z, const acb_t x, const arb_t y, long prec)

.. function:: void acb_sub(acb_t z, const acb_t x, const acb_t y, long prec)

    Sets *z* to the difference of *x* and *y*.

.. function:: void acb_mul_onei(acb_t z, const acb_t x)

    Sets *z* to *x* multiplied by the imaginary unit.

.. function:: void acb_mul_ui(acb_t z, const acb_t x, ulong y, long prec)

.. function:: void acb_mul_si(acb_t z, const acb_t x, long y, long prec)

.. function:: void acb_mul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)

.. function:: void acb_mul_arb(acb_t z, const acb_t x, const arb_t y, long prec)

    Sets *z* to the product of *x* and *y*.

.. function:: void acb_mul(acb_t z, const acb_t x, const acb_t y, long prec)

    Sets *z* to the product of *x* and *y*. If at least one part of
    *x* or *y* is zero, the operations is reduced to two real multiplications.
    If *x* and *y* are the same pointers, they are assumed to represent
    the same mathematical quantity and the squaring formula is used.

.. function:: void acb_mul_2exp_si(acb_t z, const acb_t x, long e)

.. function:: void acb_mul_2exp_fmpz(acb_t z, const acb_t x, const fmpz_t e)

    Sets *z* to *x* multiplied by `2^e`, without rounding.

.. function:: void acb_cube(acb_t z, const acb_t x, long prec)

    Sets *z* to *x* cubed, computed efficiently using two real squarings,
    two real multiplications, and scalar operations.

.. function:: void acb_addmul(acb_t z, const acb_t x, const acb_t y, long prec)

.. function:: void acb_addmul_ui(acb_t z, const acb_t x, ulong y, long prec)

.. function:: void acb_addmul_si(acb_t z, const acb_t x, long y, long prec)

.. function:: void acb_addmul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)

.. function:: void acb_addmul_arb(acb_t z, const acb_t x, const arb_t y, long prec)

    Sets *z* to *z* plus the product of *x* and *y*.

.. function:: void acb_submul(acb_t z, const acb_t x, const acb_t y, long prec)

.. function:: void acb_submul_ui(acb_t z, const acb_t x, ulong y, long prec)

.. function:: void acb_submul_si(acb_t z, const acb_t x, long y, long prec)

.. function:: void acb_submul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)

.. function:: void acb_submul_arb(acb_t z, const acb_t x, const arb_t y, long prec)

    Sets *z* to *z* minus the product of *x* and *y*.

.. function:: void acb_inv(acb_t z, const acb_t x, long prec)

    Sets *z* to the multiplicative inverse of *x*.

.. function:: void acb_div_ui(acb_t z, const acb_t x, ulong y, long prec)

.. function:: void acb_div_si(acb_t z, const acb_t x, long y, long prec)

.. function:: void acb_div_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)

.. function:: void acb_div(acb_t z, const acb_t x, const acb_t y, long prec)

    Sets *z* to the quotient of *x* and *y*.

Elementary functions
-------------------------------------------------------------------------------

.. function:: void acb_const_pi(acb_t y, long prec)

    Sets *y* to the constant `\pi`.

.. function:: void acb_log(acb_t y, const acb_t z, long prec)

    Sets *y* to the principal branch of the natural logarithm of *z*,
    computed as
    `\log(a+bi) = \frac{1}{2} \log(a^2 + b^2) + i \operatorname{arg}(a+bi)`.

.. function:: void acb_log1p(acb_t z, const acb_t x, long prec)

    Sets `z = \log(1+x)`, computed accurately when `x \approx 0`.

.. function:: void acb_exp(acb_t y, const acb_t z, long prec)

    Sets *y* to the exponential function of *z*, computed as
    `\exp(a+bi) = \exp(a) \left( \cos(b) + \sin(b) i \right)`.

.. function:: void acb_exp_pi_i(acb_t y, const acb_t z, long prec)

    Sets *y* to `\exp(\pi i z)`.

.. function:: void acb_sin(acb_t s, const acb_t z, long prec)

.. function:: void acb_cos(acb_t c, const acb_t z, long prec)

.. function:: void acb_sin_cos(arb_t s, arb_t c, const arb_t z, long prec)

    Sets `s = \sin(z)`, `c = \cos(z)`, evaluated as
    `\sin(a+bi) = \sin(a)\cosh(b) + i \cos(a)\sinh(b)`,
    `\cos(a+bi) = \cos(a)\cosh(b) - i \sin(a)\sinh(b)`.

.. function:: void acb_tan(acb_t s, const acb_t z, long prec)

    Sets `s = \tan(z) = \sin(z) / \cos(z)`, evaluated as
    `\tan(a+bi) = \sin(2a)/(\cos(2a) + \cosh(2b)) + i\sinh(2b)/(\cos(2a) + \cosh(2b))`.
    If `|b|` is small, the formula is evaluated as written; otherwise,
    we rewrite the hyperbolic functions in terms of decaying exponentials
    and evaluate the expression accurately using :func:`arb_expm1`.

.. function:: void acb_cot(acb_t s, const acb_t z, long prec)

    Sets `s = \cot(z) = \cos(z) / \sin(z)`, evaluated as
    `\cot(a+bi) = -\sin(2a)/(\cos(2a) - \cosh(2b)) + i\sinh(2b)/(\cos(2a) - \cosh(2b))`
    using the same strategy as :func:`acb_tan`.
    If `|z|` is close to zero, however, we evaluate
    `1 / \tan(z)` to avoid catastrophic cancellation.

.. function:: void acb_sin_pi(acb_t s, const acb_t z, long prec)

.. function:: void acb_cos_pi(acb_t s, const acb_t z, long prec)

.. function:: void acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, long prec)

    Sets `s = \sin(\pi z)`, `c = \cos(\pi z)`, evaluating the trigonometric
    factors of the real and imaginary part accurately via :func:`arb_sin_cos_pi`.

.. function:: void acb_tan_pi(acb_t s, const acb_t z, long prec)

    Sets `s = \tan(\pi z)`. Uses the same algorithm as :func:`acb_tan`,
    but evaluating the sine and cosine accurately via :func:`arb_sin_cos_pi`.

.. function:: void acb_cot_pi(acb_t s, const acb_t z, long prec)

    Sets `s = \cot(\pi z)`. Uses the same algorithm as :func:`acb_cot`,
    but evaluating the sine and cosine accurately via :func:`arb_sin_cos_pi`.

.. function:: void acb_atan(acb_t s, const acb_t z, long prec)

    Sets `s = \operatorname{atan}(z) = \tfrac{1}{2} i (\log(1-iz)-\log(1+iz))`.

.. function:: void acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, long prec)

.. function:: void acb_pow_ui(acb_t y, const acb_t b, ulong e, long prec)

.. function:: void acb_pow_si(acb_t y, const acb_t b, long e, long prec)

    Sets `y = b^e` using binary exponentiation (with an initial division
    if `e < 0`). Note that these functions can get slow if the exponent is
    extremely large (in such cases :func:`acb_pow` may be superior).

.. function:: void acb_pow_arb(acb_t z, const acb_t x, const arb_t y, long prec)

.. function:: void acb_pow(acb_t z, const acb_t x, const acb_t y, long prec)

    Sets `z = x^y`, computed using binary exponentiation if `y` if
    a small exact integer, as `z = (x^{1/2})^{2y}` if `y` is a small exact
    half-integer, and generally as `z = \exp(y \log x)`.

.. function:: void acb_sqrt(acb_t r, const acb_t z, long prec)

    Sets *r* to the square root of *z*.
    If either the real or imaginary part is exactly zero, only
    a single real square root is needed. Generally, we use the formula
    `\sqrt{a+bi} = u/2 + ib/u, u = \sqrt{2(|a+bi|+a)}`,
    requiring two real square root extractions.

.. function:: void acb_rsqrt(acb_t r, const acb_t z, long prec)

    Sets *r* to the reciprocal square root of *z*.
    If either the real or imaginary part is exactly zero, only
    a single real reciprocal square root is needed. Generally, we use the
    formula `1/\sqrt{a+bi} = ((a+r) - bi)/v, r = |a+bi|, v = \sqrt{r |a+bi+r|^2}`,
    requiring one real square root and one real reciprocal square root.

Rising factorials
-------------------------------------------------------------------------------

.. function:: void acb_rising_ui_bs(acb_t z, const acb_t x, ulong n, long prec)

.. function:: void acb_rising_ui_rs(acb_t z, const acb_t x, ulong n, ulong step, long prec)

.. function:: void acb_rising_ui_rec(acb_t z, const acb_t x, ulong n, long prec)

.. function:: void acb_rising_ui(acb_t z, const acb_t x, ulong n, long prec)

    Computes the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)`.

    The *bs* version uses binary splitting. The *rs* version uses rectangular
    splitting. The *rec* version uses either *bs* or *rs* depending
    on the input.
    The default version is currently identical to the *rec* version.
    In a future version, it will use the gamma function or asymptotic
    series when this is more efficient.

    The *rs* version takes an optional *step* parameter for tuning
    purposes (to use the default step length, pass zero).

.. function :: void acb_rising2_ui_bs(acb_t u, acb_t v, const acb_t x, ulong n, long prec)

.. function :: void acb_rising2_ui_rs(acb_t u, acb_t v, const acb_t x, ulong n, ulong step, long prec)

.. function :: void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, long prec)

    Letting `u(x) = x (x+1) (x+2) \cdots (x+n-1)`, simultaneously compute
    `u(x)` and `v(x) = u'(x)`, respectively using binary splitting,
    rectangular splitting (with optional nonzero step length *step*
    to override the default choice), and an automatic algorithm choice.

.. function :: void acb_rising_ui_get_mag(mag_t bound, const acb_t x, ulong n)

    Computes an upper bound for the absolute value of
    the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)`.
    Not currently optimized for large *n*.

Gamma function
-------------------------------------------------------------------------------

.. function:: void acb_gamma(acb_t y, const acb_t x, long prec)

    Computes the gamma function `y = \Gamma(x)`.

.. function:: void acb_rgamma(acb_t y, const acb_t x, long prec)

    Computes the reciprocal gamma function  `y = 1/\Gamma(x)`,
    avoiding division by zero at the poles of the gamma function.

.. function:: void acb_lgamma(acb_t y, const acb_t x, long prec)

    Computes the logarithmic gamma function `y = \log \Gamma(x)`.

    The branch cut of the logarithmic gamma function is placed on the
    negative half-axis, which means that
    `\log \Gamma(z) + \log z = \log \Gamma(z+1)` holds for all `z`,
    whereas `\log \Gamma(z) \ne \log(\Gamma(z))` in general.
    Warning: this function does not currently use the reflection
    formula, and gets very slow for `z` far into the left half-plane.

.. function:: void acb_digamma(acb_t y, const acb_t x, long prec)

    Computes the digamma function `y = \psi(x) = (\log \Gamma(x))' = \Gamma'(x) / \Gamma(x)`.

Zeta function
-------------------------------------------------------------------------------

.. function:: void acb_zeta(acb_t z, const acb_t s, long prec)

    Sets *z* to the value of the Riemann zeta function `\zeta(s)`.
    Note: for computing derivatives with respect to `s`,
    use :func:`acb_poly_zeta_series` or related methods.

.. function:: void acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, long prec)

    Sets *z* to the value of the Hurwitz zeta function `\zeta(s, a)`.
    Note: for computing derivatives with respect to `s`,
    use :func:`acb_poly_zeta_series` or related methods.

Polylogarithms
-------------------------------------------------------------------------------

.. function:: void acb_polylog(acb_t w, const acb_t s, const acb_t z, long prec)

.. function:: void acb_polylog_si(acb_t w, long s, const acb_t z, long prec)

    Sets *w* to the polylogarithm `\operatorname{Li}_s(z)`.

Arithmetic-geometric mean
-------------------------------------------------------------------------------

.. function:: void acb_agm1(acb_t m, const acb_t z, long prec)

    Sets *m* to the arithmetic-geometric mean `M(z) = \operatorname{agm}(1,z)`,
    defined such that the function is continuous in the complex plane except for
    a branch cut along the negative half axis (where it is continuous
    from above). This corresponds to always choosing an "optimal" branch for
    the square root in the arithmetic-geometric mean iteration.

.. function:: void acb_agm1_cpx(acb_ptr m, const acb_t z, long len, long prec)

    Sets the coefficients in the array *m* to the power series expansion of the
    arithmetic-geometric mean at the point *z* truncated to length *len*, i.e.
    `M(z+x) \in \mathbb{C}[[x]]`.

