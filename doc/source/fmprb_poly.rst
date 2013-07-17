**fmprb_poly.h** -- polynomials over the real numbers
===============================================================================

An :type:`fmprb_poly_t` represents a polynomial over the real numbers,
implemented as an array of coefficients of type :type:`fmprb_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (such as requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmprb_poly_struct

.. type:: fmprb_poly_t

    Contains a pointer to an array of coefficients (coeffs), the used
    length (length), and the allocated size of the array (alloc).

    An *fmprb_poly_t* is defined as an array of length one of type
    *fmprb_poly_struct*, permitting an *fmprb_poly_t* to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void fmprb_poly_init(fmprb_poly_t poly)

    Initializes the polynomial for use, setting it to the zero polynomial.

.. function:: void fmprb_poly_clear(fmprb_poly_t poly)

    Clears the polynomial, deallocating all coefficients and the
    coefficient array.

.. function:: void fmprb_poly_fit_length(fmprb_poly_t poly, long len)

    Makes sures that the coefficient array of the polynomial contains at
    least *len* initialized coefficients.

.. function:: void _fmprb_poly_set_length(fmprb_poly_t poly, long len)

    Directly changes the length of the polynomial, without allocating or
    deallocating coefficients. The value shold not exceed the allocation length.

.. function:: void _fmprb_poly_normalise(fmprb_poly_t poly)

    Strips any trailing coefficients which are identical to zero.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void fmprb_poly_zero(fmprb_poly_t poly)

.. function:: void fmprb_poly_one(fmprb_poly_t poly)

    Sets *poly* to the constant 0 respectively 1.

.. function:: void fmprb_poly_set_coeff_si(fmprb_poly_t poly, long n, long c)

.. function:: void fmprb_poly_set_coeff_fmprb(fmprb_poly_t poly, long n, const fmprb_t c)

    Sets the coefficient with index *n* in *poly* to the value *c*.
    We require that *n* is nonnegative.

.. function:: void fmprb_poly_get_coeff_fmprb(fmprb_t v, const fmprb_poly_t poly, long n)

    Sets *v* to the value of the coefficient with index *n* in *poly*.
    We require that *n* is nonnegative.

.. macro:: fmprb_poly_get_coeff_ptr(poly, n)

    Given `n >= 0`, returns a pointer to coefficient *n* of *poly*,
    or *NULL* if *n* exceeds the length of *poly*.

.. function:: void _fmprb_poly_shift_right(fmprb_struct * res, const fmprb_struct * poly, long len, long n)

.. function:: void fmprb_poly_shift_right(fmprb_poly_t res, const fmprb_poly_t poly, long n)

    Sets *res* to *poly* divided by `x^n`, throwing away the lower coefficients.
    We require that *n* is nonnegative.

.. function:: void _fmprb_poly_shift_left(fmprb_struct * res, const fmprb_struct * poly, long len, long n)

.. function:: void fmprb_poly_shift_left(fmprb_poly_t res, const fmprb_poly_t poly, long n)

    Sets *res* to *poly* multiplied by `x^n`.
    We require that *n* is nonnegative.

.. function:: void fmprb_poly_truncate(fmprb_poly_t poly, long n)

    Truncates *poly* to have length at most *n*, i.e. degree
    strictly smaller than *n*.

.. function:: long fmprb_poly_length(const fmprb_poly_t poly)

    Returns the length of *poly*, i.e. zero if *poly* is
    identically zero, and otherwise one more than the index
    of the highest term that is not identically zero.

.. function:: long fmprb_poly_degree(const fmprb_poly_t poly)

    Returns the degree of *poly*, defined as one less than its length.
    Note that if one or several leading coefficients are balls
    containing zero, this value can be larger than the true
    degree of the exact polynomial represented by *poly*,
    so the return value of this function is effectively
    an upper bound.

Conversions
-------------------------------------------------------------------------------

.. function:: void fmprb_poly_set_fmpz_poly(fmprb_poly_t poly, const fmpz_poly_t src, long prec)

.. function:: void fmprb_poly_set_fmpq_poly(fmprb_poly_t poly, const fmpq_poly_t src, long prec)

.. function:: void fmprb_poly_set_si(fmprb_poly_t poly, long src)

    Sets *poly* to *src*, rounding the coefficients to *prec* bits.


Input and output
-------------------------------------------------------------------------------

.. function:: void fmprb_poly_printd(const fmprb_poly_t poly, long digits)

    Prints the polynomial as an array of coefficients, printing each
    coefficient using *fmprb_printd*.


Random generation
-------------------------------------------------------------------------------

.. function:: void fmprb_poly_randtest(fmprb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)

    Creates a random polynomial with length at most *len*.


Comparisons
-------------------------------------------------------------------------------

.. function:: int fmprb_poly_contains_fmpq_poly(const fmprb_poly_t poly1, const fmpq_poly_t poly2)

    Returns nonzero iff *poly1* contains *poly2*.

.. function:: int fmprb_poly_equal(const fmprb_t A, const fmprb_t B)

    Returns nonzero iff *A* and *B* are equal as polynomial balls, i.e. all
    coefficients have equal midpoint and radius.

.. function:: int _fmprb_poly_overlaps(const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2)

.. function:: int fmprb_poly_overlaps(const fmprb_poly_t poly1, const fmprb_poly_t poly2)

    Returns nonzero iff *poly1* overlaps with *poly2*. The underscore
    function requires that *len1* ist at least as large as *len2*.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_add(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the sum of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void fmprb_poly_add(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long prec)

    Sets *C* to the sum of *A* and *B*.

.. function:: void _fmprb_poly_sub(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the difference of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void fmprb_poly_sub(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long prec)

    Sets *C* to the difference of *A* and *B*.

.. function:: void fmprb_poly_neg(fmprb_poly_t C, const fmprb_poly_t A)

    Sets *C* to the negation of *A*.

.. function:: void fmprb_poly_scalar_mul_2exp_si(fmprb_poly_t C, const fmprb_poly_t A, long c)

    Sets *C* to *A* multiplied by `2^c`.

.. function:: void _fmprb_poly_mullow_classical(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long n, long prec)

.. function:: void _fmprb_poly_mullow_ztrunc(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long n, long prec)

.. function:: void _fmprb_poly_mullow_block(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long n, long prec)

.. function:: void _fmprb_poly_mullow(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long n, long prec)

    Sets *{C, n}* to the product of *{A, lenA}* and *{B, lenB}*, truncated to
    length *n*. The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`,
    `n > 0`, `\mathrm{lenA} + \mathrm{lenB} - 1 \ge n`.

    The *classical* version uses a plain loop. This has good numerical
    stability but gets slow for large *n*.

    The *ztrunc* version puts each input polynomial on
    a common exponent, truncates to *prec* bits, and multiplies exactly over
    the integers. The output error is computed by cross-multiplying the
    max norms. This is fast but has poor numerical stability unless all
    coefficients are of the same magnitude.

    The *block* version decomposes the product into several
    subproducts which are computed exactly over the integers.
    This is typically nearly as fast as *ztrunc*, and the numerical
    stability is essentially as good as *classical*.

.. function:: void fmprb_poly_mullow_classical(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec)

.. function:: void fmprb_poly_mullow_ztrunc(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec)

.. function:: void fmprb_poly_mullow_block(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec)

.. function:: void fmprb_poly_mullow(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec)

    Sets *C* to the product of *A* and *B*, truncated to length *n*.

.. function:: void _fmprb_poly_mul(fmprb_struct * C, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long prec)

    Sets *{C, lenA + lenB - 1}* to the product of *{A, lenA}* and *{B, lenB}*.
    The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`.
    This function is implemented as a simple wrapper for :func:`_fmprb_poly_mullow`.

.. function:: void fmprb_poly_mul(fmprb_poly_t C, const fmprb_poly_t A, const fmprb_poly_t B, long prec)

    Sets *C* to the product of *A* and *B*.

.. function:: void _fmprb_poly_inv_series(fmprb_struct * Q, const fmprb_struct * A, long Alen, long len, long prec)

    Sets *{Q, len}* to the power series inverse of *{A, Alen}*. Uses Newton iteration.

.. function:: void fmprb_poly_inv_series(fmprb_poly_t Q, const fmprb_poly_t A, long n, long prec)

    Sets *Q* to the power series inverse of *A*, truncated to length *n*.

.. function:: void  _fmprb_poly_div_series(fmprb_struct * Q, const fmprb_struct * A, long Alen, const fmprb_struct * B, long Blen, long n, long prec)

    Sets *{Q, n}* to the power series quotient of *{A, Alen}* by *{B, Blen}*.
    Uses Newton iteration followed by multiplication.

.. function:: void fmprb_poly_div_series(fmprb_poly_t Q, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec)

    Sets *Q* to the power series quotient *A* divided by *B*, truncated to length *n*.

.. function:: void _fmprb_poly_div(fmprb_struct * Q, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long prec)

.. function:: void _fmprb_poly_rem(fmprb_struct * R, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long prec)

.. function:: void _fmprb_poly_divrem(fmprb_struct * Q, fmprb_struct * R, const fmprb_struct * A, long lenA, const fmprb_struct * B, long lenB, long prec)

.. function:: void fmprb_poly_divrem(fmprb_poly_t Q, fmprb_poly_t R, const fmprb_poly_t A, const fmprb_poly_t B, long prec)

    Performs polynomial division with remainder, computing a quotient `Q` and
    a remainder `R` such that `A = BQ + R`. The leading coefficient of `B` must
    not contain zero. The implementation reverses the inputs and performs
    power series division.

.. function:: void _fmprb_poly_div_root(fmprb_struct * Q, fmprb_t R, const fmprb_struct * A, long len, const fmprb_t c, long prec)

    Divides `A` by the polynomial `x - c`, computing the quotient `Q` as well
    as the remainder `R = f(c)`.


Composition
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_compose_horner(fmprb_struct * res, const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2, long prec)

.. function:: void fmprb_poly_compose_horner(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec)

.. function:: void _fmprb_poly_compose_divconquer(fmprb_struct * res, const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2, long prec)

.. function:: void fmprb_poly_compose_divconquer(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec)

.. function:: void _fmprb_poly_compose(fmprb_struct * res, const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2, long prec)

.. function:: void fmprb_poly_compose(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec)

    Sets *res* to the composition `h(x) = f(g(x))` where `f` is given by
    *poly1* and `g` is given by *poly2*, respectively using Horner's rule,
    divide-and-conquer, and an automatic choice between the two algorithms.
    The underscore methods do not support aliasing of the output
    with either input polynomial.

.. function:: void _fmprb_poly_compose_series_horner(fmprb_struct * res, const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2, long n, long prec)

.. function:: void fmprb_poly_compose_series_horner(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long n, long prec)

.. function:: void _fmprb_poly_compose_series_brent_kung(fmprb_struct * res, const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2, long n, long prec)

.. function:: void fmprb_poly_compose_series_brent_kung(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long n, long prec)

.. function:: void _fmprb_poly_compose_series(fmprb_struct * res, const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2, long n, long prec)

.. function:: void fmprb_poly_compose_series(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long n, long prec)

    Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*,
    respectively using Horner's rule, the Brent-Kung baby step-giant step
    algorithm, and an automatic choice between the two algorithms.
    We require that the constant term in `g(x)` is exactly zero.
    The underscore methods do not support aliasing of the output
    with either input polynomial.


.. function:: void _fmprb_poly_revert_series_lagrange(fmprb_struct * h, const fmprb_struct * f, long n, long prec)

.. function:: void fmprb_poly_revert_series_lagrange(fmprb_poly_t h, const fmprb_poly_t f, long n, long prec)

.. function:: void _fmprb_poly_revert_series_newton(fmprb_struct * h, const fmprb_struct * f, long n, long prec)

.. function:: void fmprb_poly_revert_series_newton(fmprb_poly_t h, const fmprb_poly_t f, long n, long prec)

.. function:: void _fmprb_poly_revert_series_lagrange_fast(fmprb_struct * h, const fmprb_struct * f, long n, long prec)

.. function:: void fmprb_poly_revert_series_lagrange_fast(fmprb_poly_t h, const fmprb_poly_t f, long n, long prec)

.. function:: void _fmprb_poly_revert_series(fmprb_struct * h, const fmprb_struct * f, long n, long prec)

.. function:: void fmprb_poly_revert_series(fmprb_poly_t h, const fmprb_poly_t f, long n, long prec)

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

.. function:: void _fmprb_poly_evaluate_horner(fmprb_t y, const fmprb_struct * f, long len, const fmprb_t x, long prec)

.. function:: void fmprb_poly_evaluate_horner(fmprb_t y, const fmprb_poly_t f, const fmprb_t x, long prec)

.. function:: void _fmprb_poly_evaluate_rectangular(fmprb_t y, const fmprb_struct * f, long len, const fmprb_t x, long prec)

.. function:: void fmprb_poly_evaluate_rectangular(fmprb_t y, const fmprb_poly_t f, const fmprb_t x, long prec)

.. function:: void _fmprb_poly_evaluate(fmprb_t y, const fmprb_struct * f, long len, const fmprb_t x, long prec)

.. function:: void fmprb_poly_evaluate(fmprb_t y, const fmprb_poly_t f, const fmprb_t x, long prec)

    Sets `y = f(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

.. function:: void _fmprb_poly_evaluate_fmpcb_horner(fmpcb_t y, const fmprb_struct * f, long len, const fmpcb_t x, long prec)

.. function:: void fmprb_poly_evaluate_fmpcb_horner(fmpcb_t y, const fmprb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmprb_poly_evaluate_fmpcb_rectangular(fmpcb_t y, const fmprb_struct * f, long len, const fmpcb_t x, long prec)

.. function:: void fmprb_poly_evaluate_fmpcb_rectangular(fmpcb_t y, const fmprb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmprb_poly_evaluate_fmpcb(fmpcb_t y, const fmprb_struct * f, long len, const fmpcb_t x, long prec)

.. function:: void fmprb_poly_evaluate_fmpcb(fmpcb_t y, const fmprb_poly_t f, const fmpcb_t x, long prec)

    Sets `y = f(x)` where `x` is a complex number, evaluating the
    polynomial respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

.. function:: void _fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec)

.. function:: void fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec)

.. function:: void _fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec)

.. function:: void fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec)

.. function:: void _fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec)

.. function:: void fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec)

    Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

    When Horner's rule is used, the only advantage of evaluating the
    function and its derivative simultaneously is that one does not have
    to generate the derivative polynomial explicitly.
    With the rectangular splitting algorithm, the powers can be reused,
    making simultaneous evaluation slightly faster.

.. function:: void _fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z, const fmprb_struct * f, long len, const fmpcb_t x, long prec)

.. function:: void fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmprb_poly_evaluate2_fmpcb_rectangular(fmpcb_t y, fmpcb_t z, const fmprb_struct * f, long len, const fmpcb_t x, long prec)

.. function:: void fmprb_poly_evaluate2_fmpcb_rectangular(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec)

.. function:: void _fmprb_poly_evaluate2_fmpcb(fmpcb_t y, fmpcb_t z, const fmprb_struct * f, long len, const fmpcb_t x, long prec)

.. function:: void fmprb_poly_evaluate2_fmpcb(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec)

    Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.


Product trees
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_product_roots(fmprb_struct * poly, const fmprb_struct * xs, long n, long prec)

.. function:: void fmprb_poly_product_roots(fmprb_poly_t poly, fmprb_struct * xs, long n, long prec)

    Generates the polynomial `(x-x_0)(x-x_1)\cdots(x-x_{n-1})`.

.. function:: fmprb_struct ** _fmprb_poly_tree_alloc(long len)

    Returns an initialized data structured capable of representing a
    remainder tree (product tree) of *len* roots.

.. function:: void _fmprb_poly_tree_free(fmprb_struct ** tree, long len)

    Deallocates a tree structure as allocated using *_fmprb_poly_tree_alloc*.

.. function:: void _fmprb_poly_tree_build(fmprb_struct ** tree, const fmprb_struct * roots, long len, long prec)

    Constructs a product tree from a given array of *len* roots. The tree
    structure must be pre-allocated to the specified length using
    :func:`_fmprb_poly_tree_alloc`.


Multipoint evaluation
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_evaluate_vec_iter(fmprb_struct * ys, const fmprb_struct * poly, long plen, const fmprb_struct * xs, long n, long prec)

.. function:: void fmprb_poly_evaluate_vec_iter(fmprb_struct * ys, const fmprb_poly_t poly, const fmprb_struct * xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, calling
    :func:`_fmprb_poly_evaluate` repeatedly.

.. function:: void _fmprb_poly_evaluate_vec_fast_precomp(fmprb_struct * vs, const fmprb_struct * poly, long plen, fmprb_struct ** tree, long len, long prec)

.. function:: void _fmprb_poly_evaluate_vec_fast(fmprb_struct * ys, const fmprb_struct * poly, long plen, const fmprb_struct * xs, long n, long prec)

.. function:: void fmprb_poly_evaluate_vec_fast(fmprb_struct * ys, const fmprb_poly_t poly, const fmprb_struct * xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, using
    fast multipoint evaluation.

Interpolation
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_interpolate_newton(fmprb_struct * poly, const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec)

.. function:: void fmprb_poly_interpolate_newton(fmprb_poly_t poly, const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation first interpolates in the
    Newton basis and then converts back to the monomial basis.

.. function:: void _fmprb_poly_interpolate_barycentric(fmprb_struct * poly, const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec)

.. function:: void fmprb_poly_interpolate_barycentric(fmprb_poly_t poly, const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation uses the barycentric
    form of Lagrange interpolation.

.. function:: void _fmprb_poly_interpolation_weights(fmprb_struct * w, fmprb_struct ** tree, long len, long prec)

.. function:: void _fmprb_poly_interpolate_fast_precomp(fmprb_struct * poly, const fmprb_struct * ys, fmprb_struct ** tree, const fmprb_struct * weights, long len, long prec)

.. function:: void _fmprb_poly_interpolate_fast(fmprb_struct * poly, const fmprb_struct * xs, const fmprb_struct * ys, long len, long prec)

.. function:: void fmprb_poly_interpolate_fast(fmprb_poly_t poly, const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values, using fast Lagrange interpolation.
    The precomp function takes a precomputed product tree over the
    *x* values and a vector of interpolation weights as additional inputs.


Differentiation
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_derivative(fmprb_struct * res, const fmprb_struct * poly, long len, long prec)

    Sets *{res, len - 1}* to the derivative of *{poly, len}*.
    Allows aliasing of the input and output.

.. function:: void fmprb_poly_derivative(fmprb_poly_t res, const fmprb_poly_t poly, long prec)

    Sets *res* to the derivative of *poly*.

.. function:: void _fmprb_poly_integral(fmprb_struct * res, const fmprb_struct * poly, long len, long prec)

    Sets *{res, len}* to the integral of *{poly, len - 1}*.
    Allows aliasing of the input and output.

.. function:: void fmprb_poly_integral(fmprb_poly_t res, const fmprb_poly_t poly, long prec)

    Sets *res* to the integral of *poly*.


Special functions
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_sqrt_series(fmprb_struct * g, const fmprb_struct * h, long hlen, long n, long prec)

.. function:: void fmprb_poly_sqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)

    Sets *g* to the power series square root of *h*, truncated to length *n*.
    Uses division-free Newton iteration for the reciprocal square root,
    followed by a multiplication.

    The underscore method does not support aliasing of the input and output
    arrays. It requires that *hlen* and *n* are greater than zero.

.. function:: void _fmprb_poly_rsqrt_series(fmprb_struct * g, const fmprb_struct * h, long hlen, long n, long prec)

.. function:: void fmprb_poly_rsqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)

    Sets *g* to the reciprocal power series square root of *h*, truncated to length *n*.
    Uses division-free Newton iteration.

    The underscore method does not support aliasing of the input and output
    arrays. It requires that *hlen* and *n* are greater than zero.

.. function:: void _fmprb_poly_log_series(fmprb_struct * res, const fmprb_struct * f, long flen, long n, long prec)

.. function:: void fmprb_poly_log_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)

    Sets *res* to the power series logarithm of *f*, truncated to length *n*.
    Uses the formula `\log(f(x)) = \int f'(x) / f(x) dx`, adding the logarithm of the
    constant term in *f* as the constant of integration.

    The underscore method supports aliasing of the input and output
    arrays. It requires that *flen* and *n* are greater than zero.

.. function:: void _fmprb_poly_atan_series(fmprb_struct * res, const fmprb_struct * f, long flen, long n, long prec)

.. function:: void fmprb_poly_atan_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)

.. function:: void _fmprb_poly_asin_series(fmprb_struct * res, const fmprb_struct * f, long flen, long n, long prec)

.. function:: void fmprb_poly_asin_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)

.. function:: void _fmprb_poly_acos_series(fmprb_struct * res, const fmprb_struct * f, long flen, long n, long prec)

.. function:: void fmprb_poly_acos_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)

    Sets *res* respectively to the power series inverse tangent,
    inverse sine and inverse cosine of *f*, truncated to length *n*.

    Uses the formulas

    .. math ::

        \tan^{-1}(f(x)) = \int f'(x) / (1+f(x)^2) dx,

        \sin^{-1}(f(x)) = \int f'(x) / (1-f(x)^2)^{1/2} dx,

        \cos^{-1}(f(x)) = -\int f'(x) / (1-f(x)^2)^{1/2} dx,

    adding the inverse
    function of the constant term in *f* as the constant of integration.

    The underscore methods supports aliasing of the input and output
    arrays. They require that *flen* and *n* are greater than zero.

.. function:: void _fmprb_poly_exp_series_basecase(fmprb_struct * f, const fmprb_struct * h, long hlen, long n, long prec)

.. function:: void fmprb_poly_exp_series_basecase(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec)

.. function:: void _fmprb_poly_exp_series(fmprb_struct * f, const fmprb_struct * h, long hlen, long n, long prec)

.. function:: void fmprb_poly_exp_series(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec)

    Sets `f` to the power series exponential of `h`, truncated to length `n`.

    The basecase version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The main implementation uses Newton iteration, starting from a small
    number of terms given by the basecase algorithm. The complexity
    is `O(M(n))`. Redundant operations in the Newton iteration are
    avoided by using the scheme described in [HZ2004]_.

    The underscore methods support aliasing and allow the input to be
    shorter than the output, but require the lengths to be nonzero.

.. function:: void _fmprb_poly_sin_cos_series_basecase(fmprb_struct * s, fmprb_struct * c, const fmprb_struct * h, long hlen, long n, long prec)

.. function:: void fmprb_poly_sin_cos_series_basecase(fmprb_poly_t s, fmprb_poly_t c, const fmprb_poly_t h, long n, long prec)

    Sets *s* and *c* to the power series sine and cosine of *h*.

    The basecase version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The underscore method supports aliasing and allows the input to be
    shorter than the output, but requires the lengths to be nonzero.

.. function:: void _fmprb_poly_tan_series(fmprb_struct * g, const fmprb_struct * h, long hlen, long len, long prec)

.. function:: void fmprb_poly_tan_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec)

    Sets *g* to the power series tangent of *h*.

    For small *n* takes the quotient of the sine and cosine as computed
    using the basecase algorithm. For large *n*, uses Newton iteration
    to invert the inverse tangent series. The complexity is `O(n)`.

    The underscore version does not support aliasing, and requires
    the lengths to be nonzero.

.. function:: void fmprb_poly_log_gamma_series(fmprb_poly_t f, long n, long prec)

    Sets `f` to the series expansion of `\log(\Gamma(1-x))`, truncated to
    length `n`.

.. function:: void _fmprb_poly_rfac_series_ui(fmprb_struct * res, const fmprb_struct * f, long flen, ulong r, long trunc, long prec)

.. function:: void fmprb_poly_rfac_series_ui(fmprb_poly_t res, const fmprb_poly_t f, ulong r, long trunc, long prec)

    Sets *res* to the rising factorial `(f) (f+1) (f+2) \cdots (f+r-1)`, truncated
    to length *trunc*. The underscore method assumes that *flen*, *r* and *trunc*
    are at least 1, and does not support aliasing. Uses binary splitting.


Root-finding
-------------------------------------------------------------------------------

.. function:: void _fmprb_poly_newton_convergence_factor(fmpr_t convergence_factor, const fmprb_struct * poly, long len, const fmprb_t convergence_interval, long prec)

    Given an interval `I` specified by *convergence_interval*, evaluates a bound
    for `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    where `f` is the polynomial defined by the coefficients *{poly, len}*.
    The bound is obtained by evaluating `f'(I)` and `f''(I)` directly.
    If `f` has large coefficients, `I` must be extremely precise in order to
    get a finite factor.

.. function:: int _fmprb_poly_newton_step(fmprb_t xnew, const fmprb_struct * poly, long len, const fmprb_t x, const fmprb_t convergence_interval, const fmpr_t convergence_factor, long prec)

    Performs a single step with Newton's method.

    The input consists of the polynomial `f` specified by the coefficients
    *{poly, len}*, an interval `x = [m-r, m+r]` known to contain a single root of `f`,
    an interval `I` (*convergence_interval*) containing `x` with an
    associated bound (*convergence_factor*) for
    `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    and a working precision *prec*.

    The Newton update consists of setting
    `x' = [m'-r', m'+r']` where `m' = m - f(m) / f'(m)`
    and `r' = C r^2`. The expression `m - f(m) / f'(m)` is evaluated
    using ball arithmetic at a working precision of *prec* bits, and the
    rounding error during this evaluation is accounted for in the output.
    We now check that `x' \in I` and `m' < m`. If both conditions are
    satisfied, we set *xnew* to `x'` and return nonzero.
    If either condition fails, we set *xnew* to `x` and return zero,
    indicating that no progress was made.

.. function:: void _fmprb_poly_newton_refine_root(fmprb_t r, const fmprb_struct * poly, long len, const fmprb_t start, const fmprb_t convergence_interval, const fmpr_t convergence_factor, long eval_extra_prec, long prec)

    Refines a precise estimate of a polynomial root to high precision
    by performing several Newton steps, using nearly optimally
    chosen doubling precision steps.

    The inputs are defined as for *_fmprb_poly_newton_step*, except for
    the precision parameters: *prec* is the target accuracy and
    *eval_extra_prec* is the estimated number of guard bits that need
    to be added to evaluate the polynomial accurately close to the root
    (typically, if the polynomial has large coefficients of alternating
    signs, this needs to be approximately the bit size of the coefficients).


