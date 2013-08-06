.. _fmpcb-poly:

**fmpcb_poly.h** -- polynomials over the complex numbers
===============================================================================

An :type:`fmpcb_poly_t` represents a polynomial over the complex numbers,
implemented as an array of coefficients of type :type:`fmpcb_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (such as requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpcb_poly_struct

.. type:: fmpcb_poly_t

    Contains a pointer to an array of coefficients (coeffs), the used
    length (length), and the allocated size of the array (alloc).

    An *fmpcb_poly_t* is defined as an array of length one of type
    *fmpcb_poly_struct*, permitting an *fmpcb_poly_t* to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void fmpcb_poly_init(fmpcb_poly_t poly)

    Initializes the polynomial for use, setting it to the zero polynomial.

.. function:: void fmpcb_poly_clear(fmpcb_poly_t poly)

    Clears the polynomial, deallocating all coefficients and the
    coefficient array.

.. function:: void fmpcb_poly_fit_length(fmpcb_poly_t poly, long len)

    Makes sures that the coefficient array of the polynomial contains at
    least *len* initialized coefficients.

.. function:: void _fmpcb_poly_set_length(fmpcb_poly_t poly, long len)

    Directly changes the length of the polynomial, without allocating or
    deallocating coefficients. The value shold not exceed the allocation length.

.. function:: void _fmpcb_poly_normalise(fmpcb_poly_t poly)

    Strips any trailing coefficients which are identical to zero.

.. function:: void fmpcb_poly_swap(fmpcb_poly_t poly1, fmpcb_poly_t poly2)

    Swaps *poly1* and *poly2* efficiently.


Basic properties and manipulation
-------------------------------------------------------------------------------

.. function:: long fmpcb_poly_length(const fmpcb_poly_t poly)

    Returns the length of *poly*, i.e. zero if *poly* is
    identically zero, and otherwise one more than the index
    of the highest term that is not identically zero.

.. function:: long fmpcb_poly_degree(const fmpcb_poly_t poly)

    Returns the degree of *poly*, defined as one less than its length.
    Note that if one or several leading coefficients are balls
    containing zero, this value can be larger than the true
    degree of the exact polynomial represented by *poly*,
    so the return value of this function is effectively
    an upper bound.

.. function:: void fmpcb_poly_zero(fmpcb_poly_t poly)

    Sets *poly* to the zero polynomial.

.. function:: void fmpcb_poly_one(fmpcb_poly_t poly)

    Sets *poly* to the constant polynomial 1.

.. function:: void fmpcb_poly_set(fmpcb_poly_t dest, const fmpcb_poly_t src)

    Sets *dest* to a copy of *src*.

.. function:: void fmpcb_poly_set_coeff_si(fmpcb_poly_t poly, long n, long c)

.. function:: void fmpcb_poly_set_coeff_fmpcb(fmpcb_poly_t poly, long n, const fmpcb_t c)

    Sets the coefficient with index *n* in *poly* to the value *c*.
    We require that *n* is nonnegative.

.. function:: void fmpcb_poly_get_coeff_fmpcb(fmpcb_t v, const fmpcb_poly_t poly, long n)

    Sets *v* to the value of the coefficient with index *n* in *poly*.
    We require that *n* is nonnegative.

.. macro:: fmpcb_poly_get_coeff_ptr(poly, n)

    Given `n \ge 0`, returns a pointer to coefficient *n* of *poly*,
    or *NULL* if *n* exceeds the length of *poly*.

.. function:: void _fmpcb_poly_shift_right(fmpcb_ptr res, fmpcb_srcptr poly, long len, long n)

.. function:: void fmpcb_poly_shift_right(fmpcb_poly_t res, const fmpcb_poly_t poly, long n)

    Sets *res* to *poly* divided by `x^n`, throwing away the lower coefficients.
    We require that *n* is nonnegative.

.. function:: void _fmpcb_poly_shift_left(fmpcb_ptr res, fmpcb_srcptr poly, long len, long n)

.. function:: void fmpcb_poly_shift_left(fmpcb_poly_t res, const fmpcb_poly_t poly, long n)

    Sets *res* to *poly* multiplied by `x^n`.
    We require that *n* is nonnegative.

.. function:: void fmpcb_poly_truncate(fmpcb_poly_t poly, long n)

    Truncates *poly* to have length at most *n*, i.e. degree
    strictly smaller than *n*.

Input and output
-------------------------------------------------------------------------------

.. function:: void fmpcb_poly_printd(const fmpcb_poly_t poly, long digits)

    Prints the polynomial as an array of coefficients, printing each
    coefficient using *fmprb_printd*.

Random generation
-------------------------------------------------------------------------------

.. function:: void fmpcb_poly_randtest(fmpcb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)

    Creates a random polynomial with length at most *len*.

Comparisons
-------------------------------------------------------------------------------

.. function:: int fmpcb_poly_equal(const fmpcb_poly_t A, const fmpcb_poly_t B)

    Returns nonzero iff *A* and *B* are identical as interval polynomials.

.. function:: int fmpcb_poly_contains(const fmpcb_poly_t poly1, const fmpcb_poly_t poly2)

.. function:: int fmpcb_poly_contains_fmpq_poly(const fmpcb_poly_t poly1, const fmpq_poly_t poly2)

    Returns nonzero iff *poly2* is contained in *poly1*.

.. function:: int _fmpcb_poly_overlaps(fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2)

.. function:: int fmpcb_poly_overlaps(const fmpcb_poly_t poly1, const fmpcb_poly_t poly2)

    Returns nonzero iff *poly1* overlaps with *poly2*. The underscore
    function requires that *len1* ist at least as large as *len2*.


Conversions
-------------------------------------------------------------------------------

.. function:: void fmpcb_poly_set_fmpz_poly(fmpcb_poly_t poly, const fmpz_poly_t re, long prec)

.. function:: void fmpcb_poly_set_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re)

.. function:: void fmpcb_poly_set2_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re, const fmprb_poly_t im)

.. function:: void fmpcb_poly_set_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, long prec)

.. function:: void fmpcb_poly_set2_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec)

    Sets *poly* to the given real part *re* plus the imaginary part *im*,
    both rounded to *prec* bits.

.. function:: void fmpcb_poly_set_fmpcb(fmpcb_poly_t poly, long src)

.. function:: void fmpcb_poly_set_si(fmpcb_poly_t poly, long src)

    Sets *poly* to *src*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_add(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the sum of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void fmpcb_poly_add(fmpcb_poly_t C, const fmpcb_poly_t A, const fmpcb_poly_t B, long prec)

    Sets *C* to the sum of *A* and *B*.

.. function:: void _fmpcb_poly_sub(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the difference of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void fmpcb_poly_sub(fmpcb_poly_t C, const fmpcb_poly_t A, const fmpcb_poly_t B, long prec)

    Sets *C* to the difference of *A* and *B*.

.. function:: void fmpcb_poly_neg(fmpcb_poly_t C, const fmpcb_poly_t A)

    Sets *C* to the negation of *A*.

.. function:: void fmpcb_poly_scalar_mul_2exp_si(fmpcb_poly_t C, const fmpcb_poly_t A, long c)

    Sets *C* to *A* multiplied by `2^c`.

.. function:: void _fmpcb_poly_mullow_classical(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long n, long prec)

.. function:: void _fmpcb_poly_mullow_transpose(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long n, long prec)

.. function:: void _fmpcb_poly_mullow_transpose_gauss(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long n, long prec)

.. function:: void _fmpcb_poly_mullow(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long n, long prec)

    Sets *{C, n}* to the product of *{A, lenA}* and *{B, lenB}*, truncated to
    length *n*. The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`,
    `n > 0`, `\mathrm{lenA} + \mathrm{lenB} - 1 \ge n`.

    The *classical* version uses a plain loop.

    The *transpose* version evaluates the product using four real polynomial
    multiplications (via :func:`_fmprb_poly_mullow`).

    The *transpose_gauss* version evaluates the product using three real
    polynomial multiplications. This is almost always faster than *transpose*,
    but has worse numerical stability when the coefficients vary
    in magnitude.

    The default function :func:`_fmpcb_poly_mullow` automatically switches
    been *classical* and *transpose* multiplication.

    If the input pointers are identical (and the lengths are the same),
    they are assumed to represent the same polynomial, and its
    square is computed.

.. function:: void fmpcb_poly_mullow_classical(fmpcb_poly_t C, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec)

.. function:: void fmpcb_poly_mullow_transpose(fmpcb_poly_t C, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec)

.. function:: void fmpcb_poly_mullow_transpose_gauss(fmpcb_poly_t C, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec)

.. function:: void fmpcb_poly_mullow(fmpcb_poly_t C, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec)

    Sets *C* to the product of *A* and *B*, truncated to length *n*.
    If the same variable is passed for *A* and *B*, sets *C* to the
    square of *A* truncated to length *n*.

.. function:: void _fmpcb_poly_mul(fmpcb_ptr C, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long prec)

    Sets *{C, lenA + lenB - 1}* to the product of *{A, lenA}* and *{B, lenB}*.
    The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`.
    This function is implemented as a simple wrapper for :func:`_fmpcb_poly_mullow`.

    If the input pointers are identical (and the lengths are the same),
    they are assumed to represent the same polynomial, and its
    square is computed.

.. function:: void fmpcb_poly_mul(fmpcb_poly_t C, const fmpcb_poly_t A1, const fmpcb_poly_t B2, long prec)

    Sets *C* to the product of *A* and *B*.
    If the same variable is passed for *A* and *B*, sets *C* to
    the square of *A*.

.. function:: void _fmpcb_poly_inv_series(fmpcb_ptr Qinv, fmpcb_srcptr Q, long Qlen, long len, long prec)

    Sets *{Qinv, len}* to the power series inverse of *{Q, Qlen}*. Uses Newton iteration.

.. function:: void fmpcb_poly_inv_series(fmpcb_poly_t Qinv, const fmpcb_poly_t Q, long n, long prec)

    Sets *Qinv* to the power series inverse of *Q*.

.. function:: void  _fmpcb_poly_div_series(fmpcb_ptr Q, fmpcb_srcptr A, long Alen, fmpcb_srcptr B, long Blen, long n, long prec)

    Sets *{Q, n}* to the power series quotient of *{A, Alen}* by *{B, Blen}*.
    Uses Newton iteration followed by multiplication.

.. function:: void fmpcb_poly_div_series(fmpcb_poly_t Q, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec)

    Sets *Q* to the power series quotient *A* divided by *B*, truncated to length *n*.

.. function:: void _fmpcb_poly_div(fmpcb_ptr Q, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long prec)

.. function:: void _fmpcb_poly_rem(fmpcb_ptr R, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long prec)

.. function:: void _fmpcb_poly_divrem(fmpcb_ptr Q, fmpcb_ptr R, fmpcb_srcptr A, long lenA, fmpcb_srcptr B, long lenB, long prec)

.. function:: void fmpcb_poly_divrem(fmpcb_poly_t Q, fmpcb_poly_t R, const fmpcb_poly_t A, const fmpcb_poly_t B, long prec)

    Performs polynomial division with remainder, computing a quotient `Q` and
    a remainder `R` such that `A = BQ + R`. The leading coefficient of `B` must
    not contain zero. The implementation reverses the inputs and performs
    power series division.

.. function:: void _fmpcb_poly_div_root(fmpcb_ptr Q, fmpcb_t R, fmpcb_srcptr A, long len, const fmpcb_t c, long prec)

    Divides `A` by the polynomial `x - c`, computing the quotient `Q` as well
    as the remainder `R = f(c)`.

Composition
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_compose_horner(fmpcb_ptr res, fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2, long prec)

.. function:: void fmpcb_poly_compose_horner(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long prec)

.. function:: void _fmpcb_poly_compose_divconquer(fmpcb_ptr res, fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2, long prec)

.. function:: void fmpcb_poly_compose_divconquer(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long prec)

.. function:: void _fmpcb_poly_compose(fmpcb_ptr res, fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2, long prec)

.. function:: void fmpcb_poly_compose(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long prec)

    Sets *res* to the composition `h(x) = f(g(x))` where `f` is given by
    *poly1* and `g` is given by *poly2*, respectively using Horner's rule,
    divide-and-conquer, and an automatic choice between the two algorithms.
    The underscore methods do not support aliasing of the output
    with either input polynomial.

.. function:: void _fmpcb_poly_compose_series_horner(fmpcb_ptr res, fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2, long n, long prec)

.. function:: void fmpcb_poly_compose_series_horner(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long n, long prec)

.. function:: void _fmpcb_poly_compose_series_brent_kung(fmpcb_ptr res, fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2, long n, long prec)

.. function:: void fmpcb_poly_compose_series_brent_kung(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long n, long prec)

.. function:: void _fmpcb_poly_compose_series(fmpcb_ptr res, fmpcb_srcptr poly1, long len1, fmpcb_srcptr poly2, long len2, long n, long prec)

.. function:: void fmpcb_poly_compose_series(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long n, long prec)

    Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*,
    respectively using Horner's rule, the Brent-Kung baby step-giant step
    algorithm, and an automatic choice between the two algorithms.
    We require that the constant term in `g(x)` is exactly zero.
    The underscore methods do not support aliasing of the output
    with either input polynomial.


.. function:: void _fmpcb_poly_revert_series_lagrange(fmpcb_ptr h, fmpcb_srcptr f, long n, long prec)

.. function:: void fmpcb_poly_revert_series_lagrange(fmpcb_poly_t h, const fmpcb_poly_t f, long n, long prec)

.. function:: void _fmpcb_poly_revert_series_newton(fmpcb_ptr h, fmpcb_srcptr f, long n, long prec)

.. function:: void fmpcb_poly_revert_series_newton(fmpcb_poly_t h, const fmpcb_poly_t f, long n, long prec)

.. function:: void _fmpcb_poly_revert_series_lagrange_fast(fmpcb_ptr h, fmpcb_srcptr f, long n, long prec)

.. function:: void fmpcb_poly_revert_series_lagrange_fast(fmpcb_poly_t h, const fmpcb_poly_t f, long n, long prec)

.. function:: void _fmpcb_poly_revert_series(fmpcb_ptr h, fmpcb_srcptr f, long n, long prec)

.. function:: void fmpcb_poly_revert_series(fmpcb_poly_t h, const fmpcb_poly_t f, long n, long prec)

    Sets `h` to the power series reversion of `f`, i.e. the expansion
    of the compositional inverse function `f^{-1}(x)`,
    truncated to order `O(x^n)`, using respectively
    Lagrange inversion, Newton iteration, fast Lagrange inversion,
    and a default algorithm choice.

    We require that the constant term in `f` is exactly zero and that the
    linear term is nonzero. The underscore methods assume that `f` is zero-padded to length `n`
    and do not support aliasing.

Evaluation
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_evaluate_horner(fmpcb_t y, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)

.. function:: void fmpcb_poly_evaluate_horner(fmpcb_t y, const fmpcb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmpcb_poly_evaluate_rectangular(fmpcb_t y, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)

.. function:: void fmpcb_poly_evaluate_rectangular(fmpcb_t y, const fmpcb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmpcb_poly_evaluate(fmpcb_t y, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)

.. function:: void fmpcb_poly_evaluate(fmpcb_t y, const fmpcb_poly_t f, const fmpcb_t x, long prec)

    Sets `y = f(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

.. function:: void _fmpcb_poly_evaluate2_horner(fmpcb_t y, fmpcb_t z, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)

.. function:: void fmpcb_poly_evaluate2_horner(fmpcb_t y, fmpcb_t z, const fmpcb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmpcb_poly_evaluate2_rectangular(fmpcb_t y, fmpcb_t z, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)

.. function:: void fmpcb_poly_evaluate2_rectangular(fmpcb_t y, fmpcb_t z, const fmpcb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmpcb_poly_evaluate2(fmpcb_t y, fmpcb_t z, fmpcb_srcptr f, long len, const fmpcb_t x, long prec)

.. function:: void fmpcb_poly_evaluate2(fmpcb_t y, fmpcb_t z, const fmpcb_poly_t f, const fmpcb_t x, long prec)

    Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

    When Horner's rule is used, the only advantage of evaluating the
    function and its derivative simultaneously is that one does not have
    to generate the derivative polynomial explicitly.
    With the rectangular splitting algorithm, the powers can be reused,
    making simultaneous evaluation slightly faster.


Product trees
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_product_roots(fmpcb_ptr poly, fmpcb_srcptr xs, long n, long prec)

.. function:: void fmpcb_poly_product_roots(fmpcb_poly_t poly, fmpcb_ptr xs, long n, long prec)

    Generates the polynomial `(x-x_0)(x-x_1)\cdots(x-x_{n-1})`.

.. function:: fmpcb_ptr * _fmpcb_poly_tree_alloc(long len)

    Returns an initialized data structured capable of representing a
    remainder tree (product tree) of *len* roots.

.. function:: void _fmpcb_poly_tree_free(fmpcb_ptr * tree, long len)

    Deallocates a tree structure as allocated using *_fmpcb_poly_tree_alloc*.

.. function:: void _fmpcb_poly_tree_build(fmpcb_ptr * tree, fmpcb_srcptr roots, long len, long prec)

    Constructs a product tree from a given array of *len* roots. The tree
    structure must be pre-allocated to the specified length using
    :func:`_fmpcb_poly_tree_alloc`.


Multipoint evaluation
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_evaluate_vec_iter(fmpcb_ptr ys, fmpcb_srcptr poly, long plen, fmpcb_srcptr xs, long n, long prec)

.. function:: void fmpcb_poly_evaluate_vec_iter(fmpcb_ptr ys, const fmpcb_poly_t poly, fmpcb_srcptr xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, calling
    :func:`_fmpcb_poly_evaluate` repeatedly.

.. function:: void _fmpcb_poly_evaluate_vec_fast_precomp(fmpcb_ptr vs, fmpcb_srcptr poly, long plen, fmpcb_ptr * tree, long len, long prec)

.. function:: void _fmpcb_poly_evaluate_vec_fast(fmpcb_ptr ys, fmpcb_srcptr poly, long plen, fmpcb_srcptr xs, long n, long prec)

.. function:: void fmpcb_poly_evaluate_vec_fast(fmpcb_ptr ys, const fmpcb_poly_t poly, fmpcb_srcptr xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, using
    fast multipoint evaluation.

Interpolation
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_interpolate_newton(fmpcb_ptr poly, fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)

.. function:: void fmpcb_poly_interpolate_newton(fmpcb_poly_t poly, fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation first interpolates in the
    Newton basis and then converts back to the monomial basis.

.. function:: void _fmpcb_poly_interpolate_barycentric(fmpcb_ptr poly, fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)

.. function:: void fmpcb_poly_interpolate_barycentric(fmpcb_poly_t poly, fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation uses the barycentric
    form of Lagrange interpolation.

.. function:: void _fmpcb_poly_interpolation_weights(fmpcb_ptr w, fmpcb_ptr * tree, long len, long prec)

.. function:: void _fmpcb_poly_interpolate_fast_precomp(fmpcb_ptr poly, fmpcb_srcptr ys, fmpcb_ptr * tree, fmpcb_srcptr weights, long len, long prec)

.. function:: void _fmpcb_poly_interpolate_fast(fmpcb_ptr poly, fmpcb_srcptr xs, fmpcb_srcptr ys, long len, long prec)

.. function:: void fmpcb_poly_interpolate_fast(fmpcb_poly_t poly, fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values, using fast Lagrange interpolation.
    The precomp function takes a precomputed product tree over the
    *x* values and a vector of interpolation weights as additional inputs.


Differentiation
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_derivative(fmpcb_ptr res, fmpcb_srcptr poly, long len, long prec)

    Sets *{res, len - 1}* to the derivative of *{poly, len}*.
    Allows aliasing of the input and output.

.. function:: void fmpcb_poly_derivative(fmpcb_poly_t res, const fmpcb_poly_t poly, long prec)

    Sets *res* to the derivative of *poly*.

.. function:: void _fmpcb_poly_integral(fmpcb_ptr res, fmpcb_srcptr poly, long len, long prec)

    Sets *{res, len}* to the integral of *{poly, len - 1}*.
    Allows aliasing of the input and output.

.. function:: void fmpcb_poly_integral(fmpcb_poly_t res, const fmpcb_poly_t poly, long prec)

    Sets *res* to the integral of *poly*.


Special functions
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_log_series(fmpcb_ptr res, fmpcb_srcptr f, long flen, long n, long prec)

.. function:: void fmpcb_poly_log_series(fmpcb_poly_t res, const fmpcb_poly_t f, long n, long prec)

    Sets *res* to the power series logarithm of *f*, truncated to length *n*.
    Uses the formula `\log(f(x)) = \int f'(x) / f(x) dx`, adding the logarithm of the
    constant term in *f* as the constant of integration.

    The underscore method supports aliasing of the input and output
    arrays. It requires that *flen* and *n* are greater than zero.

.. function:: void _fmpcb_poly_atan_series(fmpcb_ptr res, fmpcb_srcptr f, long flen, long n, long prec)

.. function:: void fmpcb_poly_atan_series(fmpcb_poly_t res, const fmpcb_poly_t f, long n, long prec)

    Sets *res* the power series inverse tangent of *f*, truncated to length *n*.

    Uses the formula

    .. math ::

        \tan^{-1}(f(x)) = \int f'(x) / (1+f(x)^2) dx,

    adding the function of the constant term in *f* as the constant of integration.

    The underscore method supports aliasing of the input and output
    arrays. It requires that *flen* and *n* are greater than zero.

.. function:: void _fmpcb_poly_exp_series_basecase(fmpcb_ptr f, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_exp_series_basecase(fmpcb_poly_t f, const fmpcb_poly_t h, long n, long prec)

.. function:: void _fmpcb_poly_exp_series(fmpcb_ptr f, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_exp_series(fmpcb_poly_t f, const fmpcb_poly_t h, long n, long prec)

    Sets `f` to the power series exponential of `h`, truncated to length `n`.

    The basecase version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The main implementation uses Newton iteration, starting from a small
    number of terms given by the basecase algorithm. The complexity
    is `O(M(n))`. Redundant operations in the Newton iteration are
    avoided by using the scheme described in [HZ2004]_.

    The underscore methods support aliasing and allow the input to be
    shorter than the output, but require the lengths to be nonzero.

.. function:: void _fmpcb_poly_sin_cos_series_basecase(fmpcb_ptr s, fmpcb_ptr c, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_sin_cos_series_basecase(fmpcb_poly_t s, fmpcb_poly_t c, const fmpcb_poly_t h, long n, long prec)

.. function:: void _fmpcb_poly_sin_cos_series_tangent(fmpcb_ptr s, fmpcb_ptr c, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_sin_cos_series_tangent(fmpcb_poly_t s, fmpcb_poly_t c, const fmpcb_poly_t h, long n, long prec)

.. function:: void _fmpcb_poly_sin_cos_series(fmpcb_ptr s, fmpcb_ptr c, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_sin_cos_series(fmpcb_poly_t s, fmpcb_poly_t c, const fmpcb_poly_t h, long n, long prec)

    Sets *s* and *c* to the power series sine and cosine of *h*, computed
    simultaneously.

    The *basecase* version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The *tangent* version uses the tangent half-angle formulas to compute
    the sine and cosine via :func:`_fmpcb_poly_tan_series`. This
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

.. function:: void _fmpcb_poly_sin_series(fmpcb_ptr s, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_sin_series(fmpcb_poly_t s, const fmpcb_poly_t h, long n, long prec)

.. function:: void _fmpcb_poly_cos_series(fmpcb_ptr c, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_cos_series(fmpcb_poly_t c, const fmpcb_poly_t h, long n, long prec)

    Respectively evaluates the power series sine or cosine. These functions
    simply wrap :func:`_fmpcb_poly_sin_cos_series`. The underscore methods
    support aliasing and require the lengths to be nonzero.

.. function:: void _fmpcb_poly_tan_series(fmpcb_ptr g, fmpcb_srcptr h, long hlen, long len, long prec)

.. function:: void fmpcb_poly_tan_series(fmpcb_poly_t g, const fmpcb_poly_t h, long n, long prec)

    Sets *g* to the power series tangent of *h*.

    For small *n* takes the quotient of the sine and cosine as computed
    using the basecase algorithm. For large *n*, uses Newton iteration
    to invert the inverse tangent series. The complexity is `O(M(n))`.

    The underscore version does not support aliasing, and requires
    the lengths to be nonzero.

.. function:: void _fmpcb_poly_gamma_series(fmpcb_ptr res, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_gamma_series(fmpcb_poly_t res, const fmpcb_poly_t h, long n, long prec)

.. function:: void _fmpcb_poly_rgamma_series(fmpcb_ptr res, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_rgamma_series(fmpcb_poly_t res, const fmpcb_poly_t h, long n, long prec)

.. function:: void _fmpcb_poly_lgamma_series(fmpcb_ptr res, fmpcb_srcptr h, long hlen, long n, long prec)

.. function:: void fmpcb_poly_lgamma_series(fmpcb_poly_t res, const fmpcb_poly_t h, long n, long prec)

    Sets *res* to the series expansion of `\Gamma(h(x))`, `1/\Gamma(h(x))`,
    or `\log \Gamma(h(x))`, truncated to length *n*.

    These functions first generate the Taylor series at the constant
    term of *h*, and then call :func:`_fmpcb_poly_compose_series`.
    The Taylor coefficients are generated using Stirling's series.

    The underscore methods support aliasing of the input and output
    arrays, and require that *hlen* and *n* are greater than zero.

.. function:: void _fmpcb_poly_rising_ui_series(fmpcb_ptr res, fmpcb_srcptr f, long flen, ulong r, long trunc, long prec)

.. function:: void fmpcb_poly_rising_ui_series(fmpcb_poly_t res, const fmpcb_poly_t f, ulong r, long trunc, long prec)

    Sets *res* to the rising factorial `(f) (f+1) (f+2) \cdots (f+r-1)`, truncated
    to length *trunc*. The underscore method assumes that *flen*, *r* and *trunc*
    are at least 1, and does not support aliasing. Uses binary splitting.

.. function:: void _fmpcb_poly_zeta_series(fmpcb_ptr res, fmpcb_srcptr s, long slen, const fmpcb_t a, int deflate, long n, long prec)

.. function:: void fmpcb_poly_zeta_series(fmpcb_poly_t res, const fmpcb_poly_t s, const fmpcb_t a, int deflate, long n, long prec)

    Sets *res* to the Hurwitz zeta function `\zeta(s,a)` where `s` a power
    series and `a` is a constant, truncated to length *n*.
    To evaluate the usual Riemann zeta function, set `a = 1`.

    If *deflate* is nonzero, evaluates `\zeta(s,a) + 1/(1-s)`, which
    is well-defined as a limit when the constant term of `s` is 1.
    In particular, expanding `\zeta(s,a) + 1/(1-s)` with `s = 1+x`
    gives the Stieltjes constants

    .. math ::

        \sum_{k=0}^{n-1} \frac{(-1)^k}{k!} \gamma_k(a) x^k`.

    If `a = 1`, this implementation uses the reflection formula if the midpoint
    of the constant term of `s` is negative.

Root-finding
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_root_inclusion(fmpcb_t r, const fmpcb_t m, fmpcb_srcptr poly, fmpcb_srcptr polyder, long len, long prec)

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

.. function:: long _fmpcb_poly_validate_roots(fmpcb_ptr roots, fmpcb_srcptr poly, long len, long prec)

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

.. function:: void _fmpcb_poly_refine_roots_durand_kerner(fmpcb_ptr roots, fmpcb_srcptr poly, long len, long prec)

    Refines the given roots simultaneously using a single iteration
    of the Durand-Kerner method. The radius of each root is set to an
    approximation of the correction, giving a rough estimate of its error (not
    a rigorous bound).

.. function:: long _fmpcb_poly_find_roots(fmpcb_ptr roots, fmpcb_srcptr poly, fmpcb_srcptr initial, long len, long maxiter, long prec)

.. function:: long fmpcb_poly_find_roots(fmpcb_ptr roots, const fmpcb_poly_t poly, fmpcb_srcptr initial, long maxiter, long prec)

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

