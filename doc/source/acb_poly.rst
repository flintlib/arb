.. _acb-poly:

**acb_poly.h** -- polynomials over the complex numbers
===============================================================================

An :type:`acb_poly_t` represents a polynomial over the complex numbers,
implemented as an array of coefficients of type :type:`acb_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (such as requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: acb_poly_struct

.. type:: acb_poly_t

    Contains a pointer to an array of coefficients (coeffs), the used
    length (length), and the allocated size of the array (alloc).

    An *acb_poly_t* is defined as an array of length one of type
    *acb_poly_struct*, permitting an *acb_poly_t* to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void acb_poly_init(acb_poly_t poly)

    Initializes the polynomial for use, setting it to the zero polynomial.

.. function:: void acb_poly_clear(acb_poly_t poly)

    Clears the polynomial, deallocating all coefficients and the
    coefficient array.

.. function:: void acb_poly_fit_length(acb_poly_t poly, long len)

    Makes sures that the coefficient array of the polynomial contains at
    least *len* initialized coefficients.

.. function:: void _acb_poly_set_length(acb_poly_t poly, long len)

    Directly changes the length of the polynomial, without allocating or
    deallocating coefficients. The value shold not exceed the allocation length.

.. function:: void _acb_poly_normalise(acb_poly_t poly)

    Strips any trailing coefficients which are identical to zero.

.. function:: void acb_poly_swap(acb_poly_t poly1, acb_poly_t poly2)

    Swaps *poly1* and *poly2* efficiently.


Basic properties and manipulation
-------------------------------------------------------------------------------

.. function:: long acb_poly_length(const acb_poly_t poly)

    Returns the length of *poly*, i.e. zero if *poly* is
    identically zero, and otherwise one more than the index
    of the highest term that is not identically zero.

.. function:: long acb_poly_degree(const acb_poly_t poly)

    Returns the degree of *poly*, defined as one less than its length.
    Note that if one or several leading coefficients are balls
    containing zero, this value can be larger than the true
    degree of the exact polynomial represented by *poly*,
    so the return value of this function is effectively
    an upper bound.

.. function:: void acb_poly_zero(acb_poly_t poly)

    Sets *poly* to the zero polynomial.

.. function:: void acb_poly_one(acb_poly_t poly)

    Sets *poly* to the constant polynomial 1.

.. function:: void acb_poly_set(acb_poly_t dest, const acb_poly_t src)

    Sets *dest* to a copy of *src*.

.. function:: void acb_poly_set_coeff_si(acb_poly_t poly, long n, long c)

.. function:: void acb_poly_set_coeff_acb(acb_poly_t poly, long n, const acb_t c)

    Sets the coefficient with index *n* in *poly* to the value *c*.
    We require that *n* is nonnegative.

.. function:: void acb_poly_get_coeff_acb(acb_t v, const acb_poly_t poly, long n)

    Sets *v* to the value of the coefficient with index *n* in *poly*.
    We require that *n* is nonnegative.

.. macro:: acb_poly_get_coeff_ptr(poly, n)

    Given `n \ge 0`, returns a pointer to coefficient *n* of *poly*,
    or *NULL* if *n* exceeds the length of *poly*.

.. function:: void _acb_poly_shift_right(acb_ptr res, acb_srcptr poly, long len, long n)

.. function:: void acb_poly_shift_right(acb_poly_t res, const acb_poly_t poly, long n)

    Sets *res* to *poly* divided by `x^n`, throwing away the lower coefficients.
    We require that *n* is nonnegative.

.. function:: void _acb_poly_shift_left(acb_ptr res, acb_srcptr poly, long len, long n)

.. function:: void acb_poly_shift_left(acb_poly_t res, const acb_poly_t poly, long n)

    Sets *res* to *poly* multiplied by `x^n`.
    We require that *n* is nonnegative.

.. function:: void acb_poly_truncate(acb_poly_t poly, long n)

    Truncates *poly* to have length at most *n*, i.e. degree
    strictly smaller than *n*.

Input and output
-------------------------------------------------------------------------------

.. function:: void acb_poly_printd(const acb_poly_t poly, long digits)

    Prints the polynomial as an array of coefficients, printing each
    coefficient using *arb_printd*.

Random generation
-------------------------------------------------------------------------------

.. function:: void acb_poly_randtest(acb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)

    Creates a random polynomial with length at most *len*.

Comparisons
-------------------------------------------------------------------------------

.. function:: int acb_poly_equal(const acb_poly_t A, const acb_poly_t B)

    Returns nonzero iff *A* and *B* are identical as interval polynomials.

.. function:: int acb_poly_contains(const acb_poly_t poly1, const acb_poly_t poly2)

.. function:: int acb_poly_contains_fmpz_poly(const acb_poly_t poly1, const fmpz_poly_t poly2)

.. function:: int acb_poly_contains_fmpq_poly(const acb_poly_t poly1, const fmpq_poly_t poly2)

    Returns nonzero iff *poly2* is contained in *poly1*.

.. function:: int _acb_poly_overlaps(acb_srcptr poly1, long len1, acb_srcptr poly2, long len2)

.. function:: int acb_poly_overlaps(const acb_poly_t poly1, const acb_poly_t poly2)

    Returns nonzero iff *poly1* overlaps with *poly2*. The underscore
    function requires that *len1* ist at least as large as *len2*.


Conversions
-------------------------------------------------------------------------------

.. function:: void acb_poly_set_fmpz_poly(acb_poly_t poly, const fmpz_poly_t re, long prec)

.. function:: void acb_poly_set_arb_poly(acb_poly_t poly, const arb_poly_t re)

.. function:: void acb_poly_set2_arb_poly(acb_poly_t poly, const arb_poly_t re, const arb_poly_t im)

.. function:: void acb_poly_set_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, long prec)

.. function:: void acb_poly_set2_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec)

    Sets *poly* to the given real part *re* plus the imaginary part *im*,
    both rounded to *prec* bits.

.. function:: void acb_poly_set_acb(acb_poly_t poly, long src)

.. function:: void acb_poly_set_si(acb_poly_t poly, long src)

    Sets *poly* to *src*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void _acb_poly_add(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the sum of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void acb_poly_add(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, long prec)

    Sets *C* to the sum of *A* and *B*.

.. function:: void _acb_poly_sub(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the difference of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void acb_poly_sub(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, long prec)

    Sets *C* to the difference of *A* and *B*.

.. function:: void acb_poly_neg(acb_poly_t C, const acb_poly_t A)

    Sets *C* to the negation of *A*.

.. function:: void acb_poly_scalar_mul_2exp_si(acb_poly_t C, const acb_poly_t A, long c)

    Sets *C* to *A* multiplied by `2^c`.

.. function:: void _acb_poly_mullow_classical(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long n, long prec)

.. function:: void _acb_poly_mullow_transpose(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long n, long prec)

.. function:: void _acb_poly_mullow_transpose_gauss(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long n, long prec)

.. function:: void _acb_poly_mullow(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long n, long prec)

    Sets *{C, n}* to the product of *{A, lenA}* and *{B, lenB}*, truncated to
    length *n*. The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`,
    `n > 0`, `\mathrm{lenA} + \mathrm{lenB} - 1 \ge n`.

    The *classical* version uses a plain loop.

    The *transpose* version evaluates the product using four real polynomial
    multiplications (via :func:`_arb_poly_mullow`).

    The *transpose_gauss* version evaluates the product using three real
    polynomial multiplications. This is almost always faster than *transpose*,
    but has worse numerical stability when the coefficients vary
    in magnitude.

    The default function :func:`_acb_poly_mullow` automatically switches
    been *classical* and *transpose* multiplication.

    If the input pointers are identical (and the lengths are the same),
    they are assumed to represent the same polynomial, and its
    square is computed.

.. function:: void acb_poly_mullow_classical(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, long n, long prec)

.. function:: void acb_poly_mullow_transpose(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, long n, long prec)

.. function:: void acb_poly_mullow_transpose_gauss(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, long n, long prec)

.. function:: void acb_poly_mullow(acb_poly_t C, const acb_poly_t A, const acb_poly_t B, long n, long prec)

    Sets *C* to the product of *A* and *B*, truncated to length *n*.
    If the same variable is passed for *A* and *B*, sets *C* to the
    square of *A* truncated to length *n*.

.. function:: void _acb_poly_mul(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)

    Sets *{C, lenA + lenB - 1}* to the product of *{A, lenA}* and *{B, lenB}*.
    The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`.
    This function is implemented as a simple wrapper for :func:`_acb_poly_mullow`.

    If the input pointers are identical (and the lengths are the same),
    they are assumed to represent the same polynomial, and its
    square is computed.

.. function:: void acb_poly_mul(acb_poly_t C, const acb_poly_t A1, const acb_poly_t B2, long prec)

    Sets *C* to the product of *A* and *B*.
    If the same variable is passed for *A* and *B*, sets *C* to
    the square of *A*.

.. function:: void _acb_poly_inv_series(acb_ptr Qinv, acb_srcptr Q, long Qlen, long len, long prec)

    Sets *{Qinv, len}* to the power series inverse of *{Q, Qlen}*. Uses Newton iteration.

.. function:: void acb_poly_inv_series(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)

    Sets *Qinv* to the power series inverse of *Q*.

.. function:: void  _acb_poly_div_series(acb_ptr Q, acb_srcptr A, long Alen, acb_srcptr B, long Blen, long n, long prec)

    Sets *{Q, n}* to the power series quotient of *{A, Alen}* by *{B, Blen}*.
    Uses Newton iteration followed by multiplication.

.. function:: void acb_poly_div_series(acb_poly_t Q, const acb_poly_t A, const acb_poly_t B, long n, long prec)

    Sets *Q* to the power series quotient *A* divided by *B*, truncated to length *n*.

.. function:: void _acb_poly_div(acb_ptr Q, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)

.. function:: void _acb_poly_rem(acb_ptr R, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)

.. function:: void _acb_poly_divrem(acb_ptr Q, acb_ptr R, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)

.. function:: void acb_poly_divrem(acb_poly_t Q, acb_poly_t R, const acb_poly_t A, const acb_poly_t B, long prec)

    Performs polynomial division with remainder, computing a quotient `Q` and
    a remainder `R` such that `A = BQ + R`. The leading coefficient of `B` must
    not contain zero. The implementation reverses the inputs and performs
    power series division.

.. function:: void _acb_poly_div_root(acb_ptr Q, acb_t R, acb_srcptr A, long len, const acb_t c, long prec)

    Divides `A` by the polynomial `x - c`, computing the quotient `Q` as well
    as the remainder `R = f(c)`.

Composition
-------------------------------------------------------------------------------

.. function:: void _acb_poly_compose_horner(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)

.. function:: void acb_poly_compose_horner(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)

.. function:: void _acb_poly_compose_divconquer(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)

.. function:: void acb_poly_compose_divconquer(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)

.. function:: void _acb_poly_compose(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)

.. function:: void acb_poly_compose(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)

    Sets *res* to the composition `h(x) = f(g(x))` where `f` is given by
    *poly1* and `g` is given by *poly2*, respectively using Horner's rule,
    divide-and-conquer, and an automatic choice between the two algorithms.
    The underscore methods do not support aliasing of the output
    with either input polynomial.

.. function:: void _acb_poly_compose_series_horner(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)

.. function:: void acb_poly_compose_series_horner(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)

.. function:: void _acb_poly_compose_series_brent_kung(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)

.. function:: void acb_poly_compose_series_brent_kung(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)

.. function:: void _acb_poly_compose_series(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)

.. function:: void acb_poly_compose_series(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)

    Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*,
    respectively using Horner's rule, the Brent-Kung baby step-giant step
    algorithm, and an automatic choice between the two algorithms.
    We require that the constant term in `g(x)` is exactly zero.
    The underscore methods do not support aliasing of the output
    with either input polynomial.


.. function:: void _acb_poly_revert_series_lagrange(acb_ptr h, acb_srcptr f, long flen, long n, long prec)

.. function:: void acb_poly_revert_series_lagrange(acb_poly_t h, const acb_poly_t f, long n, long prec)

.. function:: void _acb_poly_revert_series_newton(acb_ptr h, acb_srcptr f, long flen, long n, long prec)

.. function:: void acb_poly_revert_series_newton(acb_poly_t h, const acb_poly_t f, long n, long prec)

.. function:: void _acb_poly_revert_series_lagrange_fast(acb_ptr h, acb_srcptr f, long flen, long n, long prec)

.. function:: void acb_poly_revert_series_lagrange_fast(acb_poly_t h, const acb_poly_t f, long n, long prec)

.. function:: void _acb_poly_revert_series(acb_ptr h, acb_srcptr f, long flen, long n, long prec)

.. function:: void acb_poly_revert_series(acb_poly_t h, const acb_poly_t f, long n, long prec)

    Sets `h` to the power series reversion of `f`, i.e. the expansion
    of the compositional inverse function `f^{-1}(x)`,
    truncated to order `O(x^n)`, using respectively
    Lagrange inversion, Newton iteration, fast Lagrange inversion,
    and a default algorithm choice.

    We require that the constant term in `f` is exactly zero and that the
    linear term is nonzero. The underscore methods assume that *flen*
    is at least 2, and do not support aliasing.

Evaluation
-------------------------------------------------------------------------------

.. function:: void _acb_poly_evaluate_horner(acb_t y, acb_srcptr f, long len, const acb_t x, long prec)

.. function:: void acb_poly_evaluate_horner(acb_t y, const acb_poly_t f, const acb_t x, long prec)

.. function:: void _acb_poly_evaluate_rectangular(acb_t y, acb_srcptr f, long len, const acb_t x, long prec)

.. function:: void acb_poly_evaluate_rectangular(acb_t y, const acb_poly_t f, const acb_t x, long prec)

.. function:: void _acb_poly_evaluate(acb_t y, acb_srcptr f, long len, const acb_t x, long prec)

.. function:: void acb_poly_evaluate(acb_t y, const acb_poly_t f, const acb_t x, long prec)

    Sets `y = f(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

.. function:: void _acb_poly_evaluate2_horner(acb_t y, acb_t z, acb_srcptr f, long len, const acb_t x, long prec)

.. function:: void acb_poly_evaluate2_horner(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, long prec)

.. function:: void _acb_poly_evaluate2_rectangular(acb_t y, acb_t z, acb_srcptr f, long len, const acb_t x, long prec)

.. function:: void acb_poly_evaluate2_rectangular(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, long prec)

.. function:: void _acb_poly_evaluate2(acb_t y, acb_t z, acb_srcptr f, long len, const acb_t x, long prec)

.. function:: void acb_poly_evaluate2(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, long prec)

    Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

    When Horner's rule is used, the only advantage of evaluating the
    function and its derivative simultaneously is that one does not have
    to generate the derivative polynomial explicitly.
    With the rectangular splitting algorithm, the powers can be reused,
    making simultaneous evaluation slightly faster.


Product trees
-------------------------------------------------------------------------------

.. function:: void _acb_poly_product_roots(acb_ptr poly, acb_srcptr xs, long n, long prec)

.. function:: void acb_poly_product_roots(acb_poly_t poly, acb_srcptr xs, long n, long prec)

    Generates the polynomial `(x-x_0)(x-x_1)\cdots(x-x_{n-1})`.

.. function:: acb_ptr * _acb_poly_tree_alloc(long len)

    Returns an initialized data structured capable of representing a
    remainder tree (product tree) of *len* roots.

.. function:: void _acb_poly_tree_free(acb_ptr * tree, long len)

    Deallocates a tree structure as allocated using *_acb_poly_tree_alloc*.

.. function:: void _acb_poly_tree_build(acb_ptr * tree, acb_srcptr roots, long len, long prec)

    Constructs a product tree from a given array of *len* roots. The tree
    structure must be pre-allocated to the specified length using
    :func:`_acb_poly_tree_alloc`.


Multipoint evaluation
-------------------------------------------------------------------------------

.. function:: void _acb_poly_evaluate_vec_iter(acb_ptr ys, acb_srcptr poly, long plen, acb_srcptr xs, long n, long prec)

.. function:: void acb_poly_evaluate_vec_iter(acb_ptr ys, const acb_poly_t poly, acb_srcptr xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, calling
    :func:`_acb_poly_evaluate` repeatedly.

.. function:: void _acb_poly_evaluate_vec_fast_precomp(acb_ptr vs, acb_srcptr poly, long plen, acb_ptr * tree, long len, long prec)

.. function:: void _acb_poly_evaluate_vec_fast(acb_ptr ys, acb_srcptr poly, long plen, acb_srcptr xs, long n, long prec)

.. function:: void acb_poly_evaluate_vec_fast(acb_ptr ys, const acb_poly_t poly, acb_srcptr xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, using
    fast multipoint evaluation.

Interpolation
-------------------------------------------------------------------------------

.. function:: void _acb_poly_interpolate_newton(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)

.. function:: void acb_poly_interpolate_newton(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation first interpolates in the
    Newton basis and then converts back to the monomial basis.

.. function:: void _acb_poly_interpolate_barycentric(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)

.. function:: void acb_poly_interpolate_barycentric(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation uses the barycentric
    form of Lagrange interpolation.

.. function:: void _acb_poly_interpolation_weights(acb_ptr w, acb_ptr * tree, long len, long prec)

.. function:: void _acb_poly_interpolate_fast_precomp(acb_ptr poly, acb_srcptr ys, acb_ptr * tree, acb_srcptr weights, long len, long prec)

.. function:: void _acb_poly_interpolate_fast(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, long len, long prec)

.. function:: void acb_poly_interpolate_fast(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values, using fast Lagrange interpolation.
    The precomp function takes a precomputed product tree over the
    *x* values and a vector of interpolation weights as additional inputs.


Differentiation
-------------------------------------------------------------------------------

.. function:: void _acb_poly_derivative(acb_ptr res, acb_srcptr poly, long len, long prec)

    Sets *{res, len - 1}* to the derivative of *{poly, len}*.
    Allows aliasing of the input and output.

.. function:: void acb_poly_derivative(acb_poly_t res, const acb_poly_t poly, long prec)

    Sets *res* to the derivative of *poly*.

.. function:: void _acb_poly_integral(acb_ptr res, acb_srcptr poly, long len, long prec)

    Sets *{res, len}* to the integral of *{poly, len - 1}*.
    Allows aliasing of the input and output.

.. function:: void acb_poly_integral(acb_poly_t res, const acb_poly_t poly, long prec)

    Sets *res* to the integral of *poly*.


Special functions
-------------------------------------------------------------------------------

.. function:: void _acb_poly_sqrt_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_sqrt_series(acb_poly_t g, const acb_poly_t h, long n, long prec)

    Sets *g* to the power series square root of *h*, truncated to length *n*.
    Uses division-free Newton iteration for the reciprocal square root,
    followed by a multiplication.

    The underscore method does not support aliasing of the input and output
    arrays. It requires that *hlen* and *n* are greater than zero.

.. function:: void _acb_poly_rsqrt_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_rsqrt_series(acb_poly_t g, const acb_poly_t h, long n, long prec)

    Sets *g* to the reciprocal power series square root of *h*, truncated to length *n*.
    Uses division-free Newton iteration.

    The underscore method does not support aliasing of the input and output
    arrays. It requires that *hlen* and *n* are greater than zero.

.. function:: void _acb_poly_log_series(acb_ptr res, acb_srcptr f, long flen, long n, long prec)

.. function:: void acb_poly_log_series(acb_poly_t res, const acb_poly_t f, long n, long prec)

    Sets *res* to the power series logarithm of *f*, truncated to length *n*.
    Uses the formula `\log(f(x)) = \int f'(x) / f(x) dx`, adding the logarithm of the
    constant term in *f* as the constant of integration.

    The underscore method supports aliasing of the input and output
    arrays. It requires that *flen* and *n* are greater than zero.

.. function:: void _acb_poly_atan_series(acb_ptr res, acb_srcptr f, long flen, long n, long prec)

.. function:: void acb_poly_atan_series(acb_poly_t res, const acb_poly_t f, long n, long prec)

    Sets *res* the power series inverse tangent of *f*, truncated to length *n*.

    Uses the formula

    .. math ::

        \tan^{-1}(f(x)) = \int f'(x) / (1+f(x)^2) dx,

    adding the function of the constant term in *f* as the constant of integration.

    The underscore method supports aliasing of the input and output
    arrays. It requires that *flen* and *n* are greater than zero.

.. function:: void _acb_poly_exp_series_basecase(acb_ptr f, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_exp_series_basecase(acb_poly_t f, const acb_poly_t h, long n, long prec)

.. function:: void _acb_poly_exp_series(acb_ptr f, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_exp_series(acb_poly_t f, const acb_poly_t h, long n, long prec)

    Sets `f` to the power series exponential of `h`, truncated to length `n`.

    The basecase version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The main implementation uses Newton iteration, starting from a small
    number of terms given by the basecase algorithm. The complexity
    is `O(M(n))`. Redundant operations in the Newton iteration are
    avoided by using the scheme described in [HZ2004]_.

    The underscore methods support aliasing and allow the input to be
    shorter than the output, but require the lengths to be nonzero.

.. function:: void _acb_poly_sin_cos_series_basecase(acb_ptr s, acb_ptr c, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_sin_cos_series_basecase(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)

.. function:: void _acb_poly_sin_cos_series_tangent(acb_ptr s, acb_ptr c, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_sin_cos_series_tangent(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)

.. function:: void _acb_poly_sin_cos_series(acb_ptr s, acb_ptr c, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_sin_cos_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)

    Sets *s* and *c* to the power series sine and cosine of *h*, computed
    simultaneously.

    The *basecase* version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The *tangent* version uses the tangent half-angle formulas to compute
    the sine and cosine via :func:`_acb_poly_tan_series`. This
    requires `O(M(n))` operations.
    When `h = h_0 + h_1` where the constant term `h_0` is nonzero,
    the evaluation is done as
    `\sin(h_0 + h_1) = \cos(h_0) \sin(h_1) + \sin(h_0) \cos(h_1)`,
    `\cos(h_0 + h_1) = \cos(h_0) \cos(h_1) - \sin(h_0) \sin(h_1)`,
    to improve accuracy and avoid dividing by zero at the poles of
    the tangent function.

    The default version automatically selects between the *basecase* and
    *tangent* algorithms depending on the input.

    The underscore methods support aliasing and require the lengths to be nonzero.

.. function:: void _acb_poly_sin_series(acb_ptr s, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_sin_series(acb_poly_t s, const acb_poly_t h, long n, long prec)

.. function:: void _acb_poly_cos_series(acb_ptr c, acb_srcptr h, long hlen, long n, long prec)

.. function:: void acb_poly_cos_series(acb_poly_t c, const acb_poly_t h, long n, long prec)

    Respectively evaluates the power series sine or cosine. These functions
    simply wrap :func:`_acb_poly_sin_cos_series`. The underscore methods
    support aliasing and require the lengths to be nonzero.

.. function:: void _acb_poly_tan_series(acb_ptr g, acb_srcptr h, long hlen, long len, long prec)

.. function:: void acb_poly_tan_series(acb_poly_t g, const acb_poly_t h, long n, long prec)

    Sets *g* to the power series tangent of *h*.

    For small *n* takes the quotient of the sine and cosine as computed
    using the basecase algorithm. For large *n*, uses Newton iteration
    to invert the inverse tangent series. The complexity is `O(M(n))`.

    The underscore version does not support aliasing, and requires
    the lengths to be nonzero.


Root-finding
-------------------------------------------------------------------------------

.. function:: void _acb_poly_root_inclusion(acb_t r, const acb_t m, acb_srcptr poly, acb_srcptr polyder, long len, long prec)

    Given any complex number `m`, and a nonconstant polynomial `f` and its
    derivative `f'`, sets *r* to a complex interval centered on `m` that is
    guaranteed to contain at least one root of `f`.
    Such an interval is obtained by taking a ball of radius `|f(m)/f'(m)| n`
    where `n` is the degree of `f`. Proof: assume that the distance
    to the nearest root exceeds `r = |f(m)/f'(m)| n`. Then

    .. math ::

        \left|\frac{f'(m)}{f(m)}\right| =
            \left|\sum_i \frac{1}{m-\zeta_i}\right|
            \le \sum_i \frac{1}{|m-\zeta_i|}
            < \frac{n}{r} = \left|\frac{f'(m)}{f(m)}\right|

    which is a contradiction (see [Kob2010]_).

.. function:: long _acb_poly_validate_roots(acb_ptr roots, acb_srcptr poly, long len, long prec)

    Given a list of approximate roots of the input polynomial, this
    function sets a rigorous bounding interval for each root, and determines
    which roots are isolated from all the other roots.
    It then rearranges the list of roots so that the isolated roots
    are at the front of the list, and returns the count of isolated roots.

    If the return value equals the degree of the polynomial, then all
    roots have been found. If the return value is smaller, all the
    remaining output intervals are guaranteed to contain roots, but
    it is possible that not all of the polynomial's roots are contained
    among them.

.. function:: void _acb_poly_refine_roots_durand_kerner(acb_ptr roots, acb_srcptr poly, long len, long prec)

    Refines the given roots simultaneously using a single iteration
    of the Durand-Kerner method. The radius of each root is set to an
    approximation of the correction, giving a rough estimate of its error (not
    a rigorous bound).

.. function:: long _acb_poly_find_roots(acb_ptr roots, acb_srcptr poly, acb_srcptr initial, long len, long maxiter, long prec)

.. function:: long acb_poly_find_roots(acb_ptr roots, const acb_poly_t poly, acb_srcptr initial, long maxiter, long prec)

    Attempts to compute all the roots of the given nonzero polynomial *poly*
    using a working precision of *prec* bits. If *n* denotes the degree of *poly*,
    the function writes *n* approximate roots with rigorous error bounds to
    the preallocated array *roots*, and returns the number of
    roots that are isolated.

    If the return value equals the degree of the polynomial, then all
    roots have been found. If the return value is smaller, all the output
    intervals are guaranteed to contain roots, but it is possible that
    not all of the polynomial's roots are contained among them.

    The roots are computed numerically by performing several steps with
    the Durand-Kerner method and terminating if the estimated accuracy of
    the roots approaches the working precision or if the number
    of steps exceeds *maxiter*, which can be set to zero in order to use
    a default value. Finally, the approximate roots are validated rigorously.

    Initial values for the iteration can be provided as the array *initial*.
    If *initial* is set to *NULL*, default values `(0.4+0.9i)^k` are used.

    The polynomial is assumed to be squarefree. If there are repeated
    roots, the iteration is likely to find them (with low numerical accuracy),
    but the error bounds will not converge as the precision increases.

