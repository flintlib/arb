**fmpcb.h** -- complex numbers
===============================================================================

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpcb_struct

.. type:: fmpcb_t

    An *fmpcb_struct* consists of a pair of *fmprb_struct*:s.
    An *fmpcb_t* is defined as an array of length one of type
    *fmpcb_struct*, permitting an *fmpcb_t* to be passed by
    reference.

.. macro:: fmpcb_realref(x)

    Macro returning a pointer to the real part of *x* as an *fmprb_t*.

.. macro:: fmprb_imagref(x)

    Macro returning a pointer to the imaginary part of *x* as an *fmprb_t*.

Memory management
-------------------------------------------------------------------------------

.. function:: void fmpcb_init(fmprb_t x)

    Initializes the variable *x* for use. Its midpoint and radius are both
    set to zero.

.. function:: void fmpcb_clear(fmpcb_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: fmpcb_struct * _fmpcb_vec_init(long n)

    Returns a pointer to an array of *n* initialized *fmpcb_struct*:s.

.. function:: void _fmpcb_vec_clear(fmpcb_struct * v, long n)

    Clears an array of *n* initialized *fmpcb_struct*:s.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: int fmpcb_is_zero(const fmpcb_t z)

    Returns nonzero iff *z* is zero.

.. function:: int fmpcb_is_one(const fmpcb_t z)

    Returns nonzero iff *z* is exactly 1.

.. function:: int fmpcb_is_exact(const fmpcb_t z)

    Returns nonzero iff *z* is exact.

.. function:: void fmpcb_zero(fmpcb_t z)

.. function:: void fmpcb_one(fmpcb_t z)

.. function:: void fmpcb_onei(fmpcb_t z)

    Sets *z* respectively to 0, 1, `i = \sqrt{-1}`.

.. function:: void fmpcb_set(fmpcb_t z, const fmpcb_t x)

.. function:: void fmpcb_set_ui(fmpcb_t z, long x)

.. function:: void fmpcb_set_si(fmpcb_t z, long x)

.. function:: void fmpcb_set_fmpz(fmpcb_t z, const fmpz_t x)

.. function:: void fmpcb_set_fmprb(fmpcb_t z, const fmprb_t c)

    Sets *z* to the value of *x*.

.. function:: void fmpcb_set_fmpq(fmpcb_t z, const fmpq_t x, long prec)

.. function:: void fmpcb_set_round(fmpcb_t z, const fmpcb_t x, long prec)

.. function:: void fmpcb_set_round_fmpz(fmpcb_t z, const fmpz_t x, long prec)

.. function:: void fmpcb_set_round_fmprb(fmpcb_t z, const fmprb_t x, long prec)

    Sets *z* to *x*, rounded to *prec* bits.

.. function:: void fmpcb_swap(fmpcb_t z, fmpcb_t x)

    Swaps *z* and *x* efficiently.


Input and output
-------------------------------------------------------------------------------

.. function:: void fmpcb_print(const fmpcb_t x)

    Prints the internal representation of *x*.

.. function:: void fmpcb_printd(const fmpcb_t z, long digits)

    Prints *x* in decimal. The printed value of the radius is not adjusted
    to compensate for the fact that the binary-to-decimal conversion
    of both the midpoint and the radius introduces additional error.


Random number generation
-------------------------------------------------------------------------------

.. function:: void fmpcb_randtest(fmpcb_t z, flint_rand_t state, long prec, long mag_bits)

    Generates a random complex number by generating separate random
    real and imaginary parts.

Precision and comparisons
-------------------------------------------------------------------------------

.. function:: int fmpcb_equal(const fmpcb_t x, const fmpcb_t y)

    Returns nonzero iff *x* and *y* are identical.

.. function:: int fmpcb_overlaps(const fmpcb_t x, const fmpcb_t y)

    Returns nonzero iff *x* and *y* have some point in common.

.. function:: void fmpcb_get_abs_ubound_fmpr(fmpr_t u, const fmpcb_t z, long prec)

    Sets *u* to an upper bound for the absolute value of *z*, computed
    using a working precision of *prec* bits.

.. function:: void fmpcb_get_abs_lbound_fmpr(fmpr_t u, const fmpcb_t z, long prec)

    Sets *u* to a lower bound for the absolute value of *z*, computed
    using a working precision of *prec* bits.

.. function:: void fmpcb_get_rad_ubound_fmpr(fmpr_t u, const fmpcb_t z, long prec)

    Sets *u* to an upper bound for the error radius of *z* (the value
    is currently not computed tightly).

.. function:: int fmpcb_contains_fmpq(const fmpcb_t x, const fmpq_t y)

.. function:: int fmpcb_contains_fmpz(const fmpcb_t x, const fmpz_t y)

.. function:: int fmpcb_contains(const fmpcb_t x, const fmpcb_t y)

    Returns nonzero iff *y* is contained in *x*.

.. function:: int fmpcb_contains_zero(const fmpcb_t x)

    Returns nonzero iff zero is contained in *x*.

.. function:: long fmpcb_bits(const fmpcb_t x)

    Returns the maximum of *fmprb_bits* applied to the real
    and imaginary parts of *x*, i.e. the minimum precision sufficient
    to represent *x* exactly.


Complex parts
-------------------------------------------------------------------------------

.. function:: void fmpcb_arg(fmprb_t r, const fmpcb_t z, long prec)

    Sets *r* to a real interval containing the complex argument (phase) of *z*.
    We define the complex argument have a discontinuity on `(-\infty,0]`, with
    the special value `\operatorname{arg}(0) = 0`, and
    `\operatorname{arg}(a+0i) = \pi` for `a < 0`. Equivalently, if
    `z = a+bi`, the argument is given by `\operatorname{atan2}(b,a)`
    (see *fmprb_atan2*).

.. function:: void fmpcb_abs(fmprb_t r, const fmpcb_t z, long prec)

    Sets *r* to the absolute value of *z*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpcb_neg(fmpcb_t z, const fmpcb_t x)

    Sets *z* to the negation of *x*.

.. function:: void fmpcb_conj(fmpcb_t z, const fmpcb_t x)

    Sets *z* to the complex conjugate of *x*.

.. function:: void fmpcb_add_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)

.. function:: void fmpcb_add_fmpz(fmpcb_t z, const fmpcb_t x, const fmpz_t y, long prec)

.. function:: void fmpcb_add_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)

.. function:: void fmpcb_add(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

    Sets *z* to the sum of *x* and *y*.

.. function:: void fmpcb_sub_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)

.. function:: void fmpcb_sub_fmpz(fmpcb_t z, const fmpcb_t x, const fmpz_t y, long prec)

.. function:: void fmpcb_sub_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)

.. function:: void fmpcb_sub(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

    Sets *z* to the difference of *x* and *y*.

.. function:: void fmpcb_mul_onei(fmpcb_t z, const fmpcb_t x)

    Sets *z* to *x* multiplied by the imaginary unit.

.. function:: void fmpcb_mul_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)

.. function:: void fmpcb_mul_fmpz(fmpcb_t z, const fmpcb_t x, const fmpz_t y, long prec)

.. function:: void fmpcb_mul_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)

    Sets *z* to the product of *x* and *y*.

.. function:: void fmpcb_mul(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

    Sets *z* to the product of *x* and *y*. If at least one part of
    *x* or *y* is zero, the operations is reduced to two real multiplications.
    If *x* and *y* are the same pointers, they are assumed to represent
    the same mathematical quantity and the squaring formula is used.

.. function:: void fmpcb_mul_alt(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

    Sets *z* to the product of *x* and *y*. If at least one part of
    *x* or *y* is zero, the operations is reduced to two real multiplications.
    Otherwise, letting `x = a + bi`, `y = c + di`, `z = e + fi`, we use
    the formula `e = ac - bd`, `f = (a+b)(c+d) - ac - bd`,
    which requires three real multiplications instead of four.

    The drawback of this algorithm is that the numerical stability is much
    worse than for the default algorithm. In particular, if one operand
    has a large error and the other a small error, the output error will
    be about twice that of the large input error, rather than about the same.

.. function:: void fmpcb_mul_2exp_si(fmpcb_t z, const fmpcb_t x, long e)

    Sets *z* to *x* multiplied by `2^e`, without rounding.

.. function:: void fmpcb_addmul(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

.. function:: void fmpcb_addmul_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)

.. function:: void fmpcb_addmul_fmpz(fmpcb_t z, const fmpcb_t x, const fmpz_t y, long prec)

.. function:: void fmpcb_addmul_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)

    Sets *z* to *z* plus the product of *x* and *y*.

.. function:: void fmpcb_submul(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

.. function:: void fmpcb_submul_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)

.. function:: void fmpcb_submul_fmpz(fmpcb_t z, const fmpcb_t x, const fmpz_t y, long prec)

.. function:: void fmpcb_submul_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)

    Sets *z* to *z* minus the product of *x* and *y*.

.. function:: void fmpcb_inv(fmpcb_t z, const fmpcb_t x, long prec)

    Sets *z* to the multiplicative inverse of *x*.

.. function:: void fmpcb_div_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)

.. function:: void fmpcb_div_si(fmpcb_t z, const fmpcb_t x, long y, long prec)

.. function:: void fmpcb_div_fmpz(fmpcb_t z, const fmpcb_t x, const fmpz_t y, long prec)

.. function:: void fmpcb_div(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)

    Sets *z* to the quotient of *x* and *y*.

Elementary functions
-------------------------------------------------------------------------------

.. function:: void fmpcb_log(fmpcb_t y, const fmpcb_t z, long prec)

    Sets *y* to the principal branch of the natural logarithm of *z*,
    computed as
    `\log(a+bi) = \frac{1}{2} \log(a^2 + b^2) + i \operatorname{arg}(a+bi)`.

.. function:: void fmpcb_exp(fmpcb_t y, const fmpcb_t z, long prec)

    Sets *y* to the exponential function of *z*, computed as
    `\exp(a+bi) = \exp(a) \left( \cos(b) + \sin(b) i \right)`.

.. function:: void fmpcb_sin(fmpcb_t s, const fmpcb_t z, long prec)

.. function:: void fmpcb_cos(fmpcb_t c, const fmpcb_t z, long prec)

.. function:: void fmpcb_sin_cos(fmprb_t s, fmprb_t c, const fmprb_t z, long prec)

    Sets `s = \sin z`, `c = \cos z`.

.. function:: void fmpcb_sin_pi(fmpcb_t s, const fmpcb_t z, long prec)

.. function:: void fmpcb_sin_cos_pi(fmpcb_t s, fmpcb_t c, const fmpcb_t z, long prec)

    Sets `s = \sin \pi z`, `c = \cos \pi z`.

.. function:: void fmpcb_pow_fmpz(fmpcb_t y, const fmpcb_t b, const fmpz_t e, long prec)

.. function:: void fmpcb_pow_ui(fmpcb_t y, const fmpcb_t b, ulong e, long prec)

    Sets *y* to *b* raised to the power *e*, computed using binary exponentiation.

.. function:: void fmpcb_pow(fmpcb_t r, const fmpcb_t x, const fmpcb_t y, long prec)

    Sets *r* to *x* raised to the power *y*, computed as `x^y = \exp(y \log x)`.

.. function:: void fmpcb_sqrt(fmpcb_t r, const fmpcb_t z, long prec)

    Sets *r* to the square root of *z*.
    If either the real or imaginary part is exactly zero, only
    a single real square root is needed. Generally, we use the formula
    `\sqrt{a+bi} = u/2 + ib/u, u = \sqrt{2(|a+bi|+a)}`,
    requiring two real square root extractions.

.. function:: void fmpcb_invroot_newton(fmpcb_t r, const fmpcb_t a, ulong m, const fmpcb_t r0, long startprec, long prec)

    Given one inverse *m*-th root *r0* (with a valid error bound) of the complex
    number *a*, lift it from precision *startprec* to *prec* using Newton
    iteration, solving `f(z) = (1/z)^m - a = 0`.

    We require that *a* is exact and that the root is isolated from the origin.
    We also assume that the initial estimate is well isolated from the
    conjugate roots (so as to avoid converging to the wrong root).
    Given an error bound `e_n` for an input term `z_n` at step `n` of the
    Newton iteration, the error of the next term is bounded by
    `e_{n+1} < |1/f'(z)| \sum_{k=2}^{\infty} (|f^{(k)}| / k!) e_n^k`.
    Replacing `k!` by `(k-2)!` gives
    `e_{n+1} < e_n^2 (m+1) (|z_n| / (|z_n| - e_n))^{-m-2} / |z_n|`.

.. function:: void fmpcb_root_exp(fmpcb_t r, const fmpcb_t a, long m, long index, long prec)

.. function:: void fmpcb_root_newton(fmpcb_t r, const fmpcb_t a, long m, long index, long prec)

.. function:: void fmpcb_root(fmpcb_t r, const fmpcb_t a, long m, long index, long prec)

    Sets `r = \exp((1/m) (\log(a) + 2 \pi i k))`. As `k`, which is given
    by *index*, goes from `0` to `m-1`, this
    expression gives all the *m*-th roots of the complex number *a*,
    starting with the principal *m*-th root (`k = 0`).
    We allow *m* to be negative.

    The *root_exp* version evaluates the exponential directly.

    The *root_newton* version uses Newton iteration, starting from an initial
    value generated by *root_exp*. It currently requires that
    *a* is exact and requires that *m* is not equal to *LONG_MIN*.

    The *root* version makes a choice between the algorithms,
    selecting *root_newton* at high precision.

Special functions
-------------------------------------------------------------------------------

.. function:: void fmpcb_rising_ui(fmpcb_t y, const fmpcb_t x, ulong n, long prec)

    Sets *y* to the rising factorial `x (x+1) (x+2) \cdots (x+n-1)`.

.. function:: void fmpcb_gamma(fmpcb_t y, const fmpcb_t x, long prec)

    Sets `y = \Gamma(x)`, the gamma function.

.. function:: void fmpcb_rgamma(fmpcb_t y, const fmpcb_t x, long prec)

    Sets  `y = 1/\Gamma(x)`, avoiding division by zero at the poles
    of the gamma function.

.. function:: void fmpcb_lgamma(fmpcb_t y, const fmpcb_t x, long prec)

    Sets `y = \log \Gamma(x)`.

    The branch cut of the logarithmic gamma function is placed on the
    negative half-axis, which means that
    `\log \Gamma(z) + \log z = \log \Gamma(z+1)` holds for all `z`,
    whereas `\log \Gamma(z) \ne \log(\Gamma(z))` in general.
    Warning: this function does not currently use the reflection
    formula, and gets very slow for `z` far into the left half-plane.

.. function:: void fmpcb_digamma(fmpcb_t y, const fmpcb_t x, long prec)

    Sets `y = \psi(x) = (\log \Gamma(x))' = \Gamma'(x) / \Gamma(x)`.

