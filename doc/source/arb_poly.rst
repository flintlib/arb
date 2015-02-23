.. _arb-poly:

**arb_poly.h** -- polynomials over the real numbers
===============================================================================

An :type:`arb_poly_t` represents a polynomial over the real numbers,
implemented as an array of coefficients of type :type:`arb_struct`.

Most functions are provided in two versions: an underscore method which
operates directly on pre-allocated arrays of coefficients and generally
has some restrictions (such as requiring the lengths to be nonzero
and not supporting aliasing of the input and output arrays),
and a non-underscore method which performs automatic memory
management and handles degenerate cases.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: arb_poly_struct

.. type:: arb_poly_t

    Contains a pointer to an array of coefficients (coeffs), the used
    length (length), and the allocated size of the array (alloc).

    An *arb_poly_t* is defined as an array of length one of type
    *arb_poly_struct*, permitting an *arb_poly_t* to
    be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void arb_poly_init(arb_poly_t poly)

    Initializes the polynomial for use, setting it to the zero polynomial.

.. function:: void arb_poly_clear(arb_poly_t poly)

    Clears the polynomial, deallocating all coefficients and the
    coefficient array.

.. function:: void arb_poly_fit_length(arb_poly_t poly, long len)

    Makes sures that the coefficient array of the polynomial contains at
    least *len* initialized coefficients.

.. function:: void _arb_poly_set_length(arb_poly_t poly, long len)

    Directly changes the length of the polynomial, without allocating or
    deallocating coefficients. The value shold not exceed the allocation length.

.. function:: void _arb_poly_normalise(arb_poly_t poly)

    Strips any trailing coefficients which are identical to zero.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void arb_poly_zero(arb_poly_t poly)

.. function:: void arb_poly_one(arb_poly_t poly)

    Sets *poly* to the constant 0 respectively 1.

.. function:: void arb_poly_set(arb_poly_t dest, const arb_poly_t src)

    Sets *dest* to a copy of *src*.

.. function:: void arb_poly_set_round(arb_poly_t dest, const arb_poly_t src, long prec)

    Sets *dest* to a copy of *src*, rounded to *prec* bits.

.. function:: void arb_poly_set_coeff_si(arb_poly_t poly, long n, long c)

.. function:: void arb_poly_set_coeff_arb(arb_poly_t poly, long n, const arb_t c)

    Sets the coefficient with index *n* in *poly* to the value *c*.
    We require that *n* is nonnegative.

.. function:: void arb_poly_get_coeff_arb(arb_t v, const arb_poly_t poly, long n)

    Sets *v* to the value of the coefficient with index *n* in *poly*.
    We require that *n* is nonnegative.

.. macro:: arb_poly_get_coeff_ptr(poly, n)

    Given `n \ge 0`, returns a pointer to coefficient *n* of *poly*,
    or *NULL* if *n* exceeds the length of *poly*.

.. function:: void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, long len, long n)

.. function:: void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, long n)

    Sets *res* to *poly* divided by `x^n`, throwing away the lower coefficients.
    We require that *n* is nonnegative.

.. function:: void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, long len, long n)

.. function:: void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, long n)

    Sets *res* to *poly* multiplied by `x^n`.
    We require that *n* is nonnegative.

.. function:: void arb_poly_truncate(arb_poly_t poly, long n)

    Truncates *poly* to have length at most *n*, i.e. degree
    strictly smaller than *n*.

.. function:: long arb_poly_length(const arb_poly_t poly)

    Returns the length of *poly*, i.e. zero if *poly* is
    identically zero, and otherwise one more than the index
    of the highest term that is not identically zero.

.. function:: long arb_poly_degree(const arb_poly_t poly)

    Returns the degree of *poly*, defined as one less than its length.
    Note that if one or several leading coefficients are balls
    containing zero, this value can be larger than the true
    degree of the exact polynomial represented by *poly*,
    so the return value of this function is effectively
    an upper bound.

Conversions
-------------------------------------------------------------------------------

.. function:: void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, long prec)

.. function:: void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, long prec)

.. function:: void arb_poly_set_si(arb_poly_t poly, long src)

    Sets *poly* to *src*, rounding the coefficients to *prec* bits.


Input and output
-------------------------------------------------------------------------------

.. function:: void arb_poly_printd(const arb_poly_t poly, long digits)

    Prints the polynomial as an array of coefficients, printing each
    coefficient using *arb_printd*.


Random generation
-------------------------------------------------------------------------------

.. function:: void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)

    Creates a random polynomial with length at most *len*.


Comparisons
-------------------------------------------------------------------------------

.. function:: int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2)

.. function:: int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2)

.. function:: int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2)

    Returns nonzero iff *poly1* contains *poly2*.

.. function:: int arb_poly_equal(const arb_poly_t A, const arb_poly_t B)

    Returns nonzero iff *A* and *B* are equal as polynomial balls, i.e. all
    coefficients have equal midpoint and radius.

.. function:: int _arb_poly_overlaps(arb_srcptr poly1, long len1, arb_srcptr poly2, long len2)

.. function:: int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2)

    Returns nonzero iff *poly1* overlaps with *poly2*. The underscore
    function requires that *len1* ist at least as large as *len2*.

.. function:: int arb_poly_get_unique_fmpz_poly(fmpz_poly_t z, const arb_poly_t x)

    If *x* contains a unique integer polynomial, sets *z* to that value and returns
    nonzero. Otherwise (if *x* represents no integers or more than one integer),
    returns zero, possibly partially modifying *z*.

Bounds
-------------------------------------------------------------------------------

.. function:: void _arb_poly_majorant(arb_ptr res, arb_srcptr poly, long len, long prec)

.. function:: void arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, long prec)

    Sets *res* to an exact real polynomial whose coefficients are
    upper bounds for the absolute values of the coefficients in *poly*,
    rounded to *prec* bits.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void _arb_poly_add(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the sum of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void arb_poly_add(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long prec)

.. function:: void arb_poly_add_si(arb_poly_t C, const arb_poly_t A, long B, long prec)

    Sets *C* to the sum of *A* and *B*.

.. function:: void _arb_poly_sub(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

    Sets *{C, max(lenA, lenB)}* to the difference of *{A, lenA}* and *{B, lenB}*.
    Allows aliasing of the input and output operands.

.. function:: void arb_poly_sub(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long prec)

    Sets *C* to the difference of *A* and *B*.

.. function:: void arb_poly_neg(arb_poly_t C, const arb_poly_t A)

    Sets *C* to the negation of *A*.

.. function:: void arb_poly_scalar_mul_2exp_si(arb_poly_t C, const arb_poly_t A, long c)

    Sets *C* to *A* multiplied by `2^c`.

.. function:: void _arb_poly_mullow_classical(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)

.. function:: void _arb_poly_mullow_block(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)

.. function:: void _arb_poly_mullow(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)

    Sets *{C, n}* to the product of *{A, lenA}* and *{B, lenB}*, truncated to
    length *n*. The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`,
    `n > 0`, `\mathrm{lenA} + \mathrm{lenB} - 1 \ge n`.

    The *classical* version uses a plain loop. This has good numerical
    stability but gets slow for large *n*.

    The *block* version decomposes the product into several
    subproducts which are computed exactly over the integers.

    It first attempts to find an integer `c`
    such that `A(2^c x)` and `B(2^c x)` have slowly varying
    coefficients, to reduce the number of blocks.

    The scaling factor `c` is chosen in a quick, heuristic way
    by picking the first and last nonzero terms in each polynomial.
    If the indices in `A` are `a_2, a_1` and the log-2 magnitudes
    are `e_2, e_1`, and the indices in `B` are `b_2, b_1`
    with corresponding magnitudes `f_2, f_1`, then we compute
    `c` as the weighted arithmetic mean of the slopes,
    rounded to the nearest integer:

    .. math ::

        c = \left\lfloor
            \frac{(e_2 - e_1) + (f_2 + f_1)}{(a_2 - a_1) + (b_2 - b_1)}
            + \frac{1}{2}
            \right \rfloor.

    This strategy is used because it is simple. It is not optimal
    in all cases, but will typically give good performance when
    multiplying two power series with a similar decay rate.

    The default algorithm chooses the *classical* algorithm for
    short polynomials and the *block* algorithm for long polynomials.

    If the input pointers are identical (and the lengths are the same),
    they are assumed to represent the same polynomial, and its
    square is computed.

.. function:: void arb_poly_mullow_classical(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

.. function:: void arb_poly_mullow_ztrunc(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

.. function:: void arb_poly_mullow_block(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

.. function:: void arb_poly_mullow(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long n, long prec)

    Sets *C* to the product of *A* and *B*, truncated to length *n*.
    If the same variable is passed for *A* and *B*, sets *C* to the square
    of *A* truncated to length *n*.

.. function:: void _arb_poly_mul(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

    Sets *{C, lenA + lenB - 1}* to the product of *{A, lenA}* and *{B, lenB}*.
    The output is not allowed to be aliased with either of the
    inputs. We require `\mathrm{lenA} \ge \mathrm{lenB} > 0`.
    This function is implemented as a simple wrapper for :func:`_arb_poly_mullow`.

    If the input pointers are identical (and the lengths are the same),
    they are assumed to represent the same polynomial, and its
    square is computed.

.. function:: void arb_poly_mul(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, long prec)

    Sets *C* to the product of *A* and *B*.
    If the same variable is passed for *A* and *B*, sets *C* to the
    square of *A*.

.. function:: void _arb_poly_inv_series(arb_ptr Q, arb_srcptr A, long Alen, long len, long prec)

    Sets *{Q, len}* to the power series inverse of *{A, Alen}*. Uses Newton iteration.

.. function:: void arb_poly_inv_series(arb_poly_t Q, const arb_poly_t A, long n, long prec)

    Sets *Q* to the power series inverse of *A*, truncated to length *n*.

.. function:: void  _arb_poly_div_series(arb_ptr Q, arb_srcptr A, long Alen, arb_srcptr B, long Blen, long n, long prec)

    Sets *{Q, n}* to the power series quotient of *{A, Alen}* by *{B, Blen}*.
    Uses Newton iteration followed by multiplication.

.. function:: void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, long n, long prec)

    Sets *Q* to the power series quotient *A* divided by *B*, truncated to length *n*.

.. function:: void _arb_poly_div(arb_ptr Q, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

.. function:: void _arb_poly_rem(arb_ptr R, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

.. function:: void _arb_poly_divrem(arb_ptr Q, arb_ptr R, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)

.. function:: int arb_poly_divrem(arb_poly_t Q, arb_poly_t R, const arb_poly_t A, const arb_poly_t B, long prec)

    Performs polynomial division with remainder, computing a quotient `Q` and
    a remainder `R` such that `A = BQ + R`. The implementation reverses the
    inputs and performs power series division.

    If the leading coefficient of `B` contains zero (or if `B` is identically
    zero), returns 0 indicating failure without modifying the outputs.
    Otherwise returns nonzero.

.. function:: void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A, long len, const arb_t c, long prec)

    Divides `A` by the polynomial `x - c`, computing the quotient `Q` as well
    as the remainder `R = f(c)`.


Composition
-------------------------------------------------------------------------------

.. function:: void _arb_poly_compose_horner(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)

.. function:: void arb_poly_compose_horner(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)

.. function:: void _arb_poly_compose_divconquer(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)

.. function:: void arb_poly_compose_divconquer(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)

.. function:: void _arb_poly_compose(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)

.. function:: void arb_poly_compose(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)

    Sets *res* to the composition `h(x) = f(g(x))` where `f` is given by
    *poly1* and `g` is given by *poly2*, respectively using Horner's rule,
    divide-and-conquer, and an automatic choice between the two algorithms.
    The underscore methods do not support aliasing of the output
    with either input polynomial.

.. function:: void _arb_poly_compose_series_horner(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)

.. function:: void arb_poly_compose_series_horner(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)

.. function:: void _arb_poly_compose_series_brent_kung(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)

.. function:: void arb_poly_compose_series_brent_kung(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)

.. function:: void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)

.. function:: void arb_poly_compose_series(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)

    Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*,
    respectively using Horner's rule, the Brent-Kung baby step-giant step
    algorithm, and an automatic choice between the two algorithms.
    We require that the constant term in `g(x)` is exactly zero.
    The underscore methods do not support aliasing of the output
    with either input polynomial.


.. function:: void _arb_poly_revert_series_lagrange(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_revert_series_lagrange(arb_poly_t h, const arb_poly_t f, long n, long prec)

.. function:: void _arb_poly_revert_series_newton(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_revert_series_newton(arb_poly_t h, const arb_poly_t f, long n, long prec)

.. function:: void _arb_poly_revert_series_lagrange_fast(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_revert_series_lagrange_fast(arb_poly_t h, const arb_poly_t f, long n, long prec)

.. function:: void _arb_poly_revert_series(arb_ptr h, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_revert_series(arb_poly_t h, const arb_poly_t f, long n, long prec)

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

.. function:: void _arb_poly_evaluate_horner(arb_t y, arb_srcptr f, long len, const arb_t x, long prec)

.. function:: void arb_poly_evaluate_horner(arb_t y, const arb_poly_t f, const arb_t x, long prec)

.. function:: void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr f, long len, const arb_t x, long prec)

.. function:: void arb_poly_evaluate_rectangular(arb_t y, const arb_poly_t f, const arb_t x, long prec)

.. function:: void _arb_poly_evaluate(arb_t y, arb_srcptr f, long len, const arb_t x, long prec)

.. function:: void arb_poly_evaluate(arb_t y, const arb_poly_t f, const arb_t x, long prec)

    Sets `y = f(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

.. function:: void _arb_poly_evaluate_acb_horner(acb_t y, arb_srcptr f, long len, const acb_t x, long prec)

.. function:: void arb_poly_evaluate_acb_horner(acb_t y, const arb_poly_t f, const acb_t x, long prec)

.. function:: void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr f, long len, const acb_t x, long prec)

.. function:: void arb_poly_evaluate_acb_rectangular(acb_t y, const arb_poly_t f, const acb_t x, long prec)

.. function:: void _arb_poly_evaluate_acb(acb_t y, arb_srcptr f, long len, const acb_t x, long prec)

.. function:: void arb_poly_evaluate_acb(acb_t y, const arb_poly_t f, const acb_t x, long prec)

    Sets `y = f(x)` where `x` is a complex number, evaluating the
    polynomial respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

.. function:: void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)

.. function:: void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)

.. function:: void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)

.. function:: void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)

.. function:: void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)

.. function:: void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)

    Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.

    When Horner's rule is used, the only advantage of evaluating the
    function and its derivative simultaneously is that one does not have
    to generate the derivative polynomial explicitly.
    With the rectangular splitting algorithm, the powers can be reused,
    making simultaneous evaluation slightly faster.

.. function:: void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)

.. function:: void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)

.. function:: void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)

.. function:: void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)

.. function:: void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)

.. function:: void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)

    Sets `y = f(x), z = f'(x)`, evaluated respectively using Horner's rule,
    rectangular splitting, and an automatic algorithm choice.


Product trees
-------------------------------------------------------------------------------

.. function:: void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, long n, long prec)

.. function:: void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, long n, long prec)

    Generates the polynomial `(x-x_0)(x-x_1)\cdots(x-x_{n-1})`.

.. function:: arb_ptr * _arb_poly_tree_alloc(long len)

    Returns an initialized data structured capable of representing a
    remainder tree (product tree) of *len* roots.

.. function:: void _arb_poly_tree_free(arb_ptr * tree, long len)

    Deallocates a tree structure as allocated using *_arb_poly_tree_alloc*.

.. function:: void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, long len, long prec)

    Constructs a product tree from a given array of *len* roots. The tree
    structure must be pre-allocated to the specified length using
    :func:`_arb_poly_tree_alloc`.


Multipoint evaluation
-------------------------------------------------------------------------------

.. function:: void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, long plen, arb_srcptr xs, long n, long prec)

.. function:: void arb_poly_evaluate_vec_iter(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, calling
    :func:`_arb_poly_evaluate` repeatedly.

.. function:: void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly, long plen, arb_ptr * tree, long len, long prec)

.. function:: void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, long plen, arb_srcptr xs, long n, long prec)

.. function:: void arb_poly_evaluate_vec_fast(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, long n, long prec)

    Evaluates the polynomial simultaneously at *n* given points, using
    fast multipoint evaluation.

Interpolation
-------------------------------------------------------------------------------

.. function:: void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

.. function:: void arb_poly_interpolate_newton(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation first interpolates in the
    Newton basis and then converts back to the monomial basis.

.. function:: void _arb_poly_interpolate_barycentric(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

.. function:: void arb_poly_interpolate_barycentric(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values. This implementation uses the barycentric
    form of Lagrange interpolation.

.. function:: void _arb_poly_interpolation_weights(arb_ptr w, arb_ptr * tree, long len, long prec)

.. function:: void _arb_poly_interpolate_fast_precomp(arb_ptr poly, arb_srcptr ys, arb_ptr * tree, arb_srcptr weights, long len, long prec)

.. function:: void _arb_poly_interpolate_fast(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long len, long prec)

.. function:: void arb_poly_interpolate_fast(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

    Recovers the unique polynomial of length at most *n* that interpolates
    the given *x* and *y* values, using fast Lagrange interpolation.
    The precomp function takes a precomputed product tree over the
    *x* values and a vector of interpolation weights as additional inputs.


Differentiation
-------------------------------------------------------------------------------

.. function:: void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, long len, long prec)

    Sets *{res, len - 1}* to the derivative of *{poly, len}*.
    Allows aliasing of the input and output.

.. function:: void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, long prec)

    Sets *res* to the derivative of *poly*.

.. function:: void _arb_poly_integral(arb_ptr res, arb_srcptr poly, long len, long prec)

    Sets *{res, len}* to the integral of *{poly, len - 1}*.
    Allows aliasing of the input and output.

.. function:: void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, long prec)

    Sets *res* to the integral of *poly*.


Transforms
-------------------------------------------------------------------------------

.. function:: void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec)

.. function:: void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec)

    Computes the Borel transform of the input polynomial, mapping `\sum_k a_k x^k`
    to `\sum_k (a_k / k!) x^k`. The underscore method allows aliasing.

.. function:: void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec)

.. function:: void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec)

    Computes the inverse Borel transform of the input polynomial, mapping `\sum_k a_k x^k`
    to `\sum_k a_k k! x^k`. The underscore method allows aliasing.

.. function:: void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, long alen, long len, long prec)

.. function:: void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, long len, long prec)

.. function:: void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, long alen, long len, long prec)

.. function:: void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, long len, long prec)

.. function:: void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, long alen, long len, long prec)

.. function:: void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, long len, long prec)

    Computes the binomial transform of the input polynomial, truncating
    the output to length *len*.
    The binomial transform maps the coefficients `a_k` in the input polynomial
    to the coefficients `b_k` in the output polynomial via
    `b_n = \sum_{k=0}^n (-1)^k {n \choose k} a_k`.
    The binomial transform is equivalent to the power series composition
    `f(x) \to (1-x)^{-1} f(x/(x-1))`, and is its own inverse.

    The *basecase* version evaluates coefficients one by one from the
    definition, generating the binomial coefficients by a recurrence
    relation.

    The *convolution* version uses the identity
    `T(f(x)) = B^{-1}(e^x B(f(-x)))` where `T` denotes the binomial
    transform operator and `B` denotes the Borel transform operator.
    This only costs a single polynomial multiplication, plus some
    scalar operations.

    The default version automatically chooses an algorithm.

    The underscore methods do not support aliasing, and assume that
    the lengths are nonzero.

Powers and elementary functions
-------------------------------------------------------------------------------

.. function:: void _arb_poly_pow_ui_trunc_binexp(arb_ptr res, arb_srcptr f, long flen, ulong exp, long len, long prec)

    Sets *{res, len}* to *{f, flen}* raised to the power *exp*, truncated
    to length *len*. Requires that *len* is no longer than the length
    of the power as computed without truncation (i.e. no zero-padding is performed).
    Does not support aliasing of the input and output, and requires
    that *flen* and *len* are positive.
    Uses binary expontiation.

.. function:: void arb_poly_pow_ui_trunc_binexp(arb_poly_t res, const arb_poly_t poly, ulong exp, long len, long prec)

    Sets *res* to *poly* raised to the power *exp*, truncated to length *len*.
    Uses binary exponentiation.

.. function:: void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, long flen, ulong exp, long prec)

    Sets *res* to *{f, flen}* raised to the power *exp*. Does not
    support aliasing of the input and output, and requires that
    *flen* is positive.

.. function:: void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, ulong exp, long prec)

    Sets *res* to *poly* raised to the power *exp*.

.. function:: void _arb_poly_pow_series(arb_ptr h, arb_srcptr f, long flen, arb_srcptr g, long glen, long len, long prec)

    Sets *{h, len}* to the power series `f(x)^{g(x)} = \exp(g(x) \log f(x))` truncated
    to length *len*. This function detects special cases such as *g* being an
    exact small integer or `\pm 1/2`, and computes such powers more
    efficiently. This function does not support aliasing of the output
    with either of the input operands. It requires that all lengths
    are positive, and assumes that *flen* and *glen* do not exceed *len*.

.. function:: void arb_poly_pow_series(arb_poly_t h, const arb_poly_t f, const arb_poly_t g, long len, long prec)

    Sets *h* to the power series `f(x)^{g(x)} = \exp(g(x) \log f(x))` truncated
    to length *len*. This function detects special cases such as *g* being an
    exact small integer or `\pm 1/2`, and computes such powers more
    efficiently.

.. function:: void _arb_poly_pow_arb_series(arb_ptr h, arb_srcptr f, long flen, const arb_t g, long len, long prec)

    Sets *{h, len}* to the power series `f(x)^g = \exp(g \log f(x))` truncated
    to length *len*. This function detects special cases such as *g* being an
    exact small integer or `\pm 1/2`, and computes such powers more
    efficiently. This function does not support aliasing of the output
    with either of the input operands. It requires that all lengths
    are positive, and assumes that *flen* does not exceed *len*.

.. function:: void arb_poly_pow_arb_series(arb_poly_t h, const arb_poly_t f, const arb_t g, long len, long prec)

    Sets *h* to the power series `f(x)^g = \exp(g \log f(x))` truncated
    to length *len*.

.. function:: void _arb_poly_sqrt_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)

    Sets *g* to the power series square root of *h*, truncated to length *n*.
    Uses division-free Newton iteration for the reciprocal square root,
    followed by a multiplication.

    The underscore method does not support aliasing of the input and output
    arrays. It requires that *hlen* and *n* are greater than zero.

.. function:: void _arb_poly_rsqrt_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)

    Sets *g* to the reciprocal power series square root of *h*, truncated to length *n*.
    Uses division-free Newton iteration.

    The underscore method does not support aliasing of the input and output
    arrays. It requires that *hlen* and *n* are greater than zero.

.. function:: void _arb_poly_log_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, long n, long prec)

    Sets *res* to the power series logarithm of *f*, truncated to length *n*.
    Uses the formula `\log(f(x)) = \int f'(x) / f(x) dx`, adding the logarithm of the
    constant term in *f* as the constant of integration.

    The underscore method supports aliasing of the input and output
    arrays. It requires that *flen* and *n* are greater than zero.

.. function:: void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, long n, long prec)

.. function:: void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, long n, long prec)

.. function:: void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)

.. function:: void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, long n, long prec)

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

.. function:: void _arb_poly_exp_series_basecase(arb_ptr f, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, long n, long prec)

.. function:: void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, long n, long prec)

    Sets `f` to the power series exponential of `h`, truncated to length `n`.

    The basecase version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The main implementation uses Newton iteration, starting from a small
    number of terms given by the basecase algorithm. The complexity
    is `O(M(n))`. Redundant operations in the Newton iteration are
    avoided by using the scheme described in [HZ2004]_.

    The underscore methods support aliasing and allow the input to be
    shorter than the output, but require the lengths to be nonzero.

.. function:: void _arb_poly_sin_cos_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec, int times_pi)

.. function:: void arb_poly_sin_cos_series_basecase(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec, int times_pi)

.. function:: void _arb_poly_sin_cos_series_tangent(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec, int times_pi)

.. function:: void arb_poly_sin_cos_series_tangent(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec, int times_pi)

.. function:: void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)

    Sets *s* and *c* to the power series sine and cosine of *h*, computed
    simultaneously.

    The *basecase* version uses a simple recurrence for the coefficients,
    requiring `O(nm)` operations where `m` is the length of `h`.

    The *tangent* version uses the tangent half-angle formulas to compute
    the sine and cosine via :func:`_arb_poly_tan_series`. This
    requires `O(M(n))` operations.
    When `h = h_0 + h_1` where the constant term `h_0` is nonzero,
    the evaluation is done as
    `\sin(h_0 + h_1) = \cos(h_0) \sin(h_1) + \sin(h_0) \cos(h_1)`,
    `\cos(h_0 + h_1) = \cos(h_0) \cos(h_1) - \sin(h_0) \sin(h_1)`,
    to improve accuracy and avoid dividing by zero at the poles of
    the tangent function.

    The default version automatically selects between the *basecase* and
    *tangent* algorithms depending on the input.

    The *basecase* and *tangent* versions take a flag *times_pi*
    specifying that the input is to be multiplied by `\pi`.

    The underscore methods support aliasing and require the lengths to be nonzero.

.. function:: void _arb_poly_sin_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_sin_series(arb_poly_t s, const arb_poly_t h, long n, long prec)

.. function:: void _arb_poly_cos_series(arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_cos_series(arb_poly_t c, const arb_poly_t h, long n, long prec)

    Respectively evaluates the power series sine or cosine. These functions
    simply wrap :func:`_arb_poly_sin_cos_series`. The underscore methods
    support aliasing and require the lengths to be nonzero.

.. function:: void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)

.. function:: void _arb_poly_sin_pi_series(arb_ptr s, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_sin_pi_series(arb_poly_t s, const arb_poly_t h, long n, long prec)

.. function:: void _arb_poly_cos_pi_series(arb_ptr c, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_cos_pi_series(arb_poly_t c, const arb_poly_t h, long n, long prec)

    Compute the respective trigonometric functions of the input
    multiplied by `\pi`.

.. function:: void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, long hlen, long len, long prec)

.. function:: void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, long n, long prec)

    Sets *g* to the power series tangent of *h*.

    For small *n* takes the quotient of the sine and cosine as computed
    using the basecase algorithm. For large *n*, uses Newton iteration
    to invert the inverse tangent series. The complexity is `O(M(n))`.

    The underscore version does not support aliasing, and requires
    the lengths to be nonzero.

Gamma function and factorials
-------------------------------------------------------------------------------

.. function:: void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

.. function:: void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

.. function:: void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    Sets *res* to the series expansion of `\Gamma(h(x))`, `1/\Gamma(h(x))`,
    or `\log \Gamma(h(x))`, truncated to length *n*.

    These functions first generate the Taylor series at the constant
    term of *h*, and then call :func:`_arb_poly_compose_series`.
    The Taylor coefficients are generated using the Riemann zeta function
    if the constant term of *h* is a small integer,
    and with Stirling's series otherwise.

    The underscore methods support aliasing of the input and output
    arrays, and require that *hlen* and *n* are greater than zero.

.. function:: void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, long flen, ulong r, long trunc, long prec)

.. function:: void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, long trunc, long prec)

    Sets *res* to the rising factorial `(f) (f+1) (f+2) \cdots (f+r-1)`, truncated
    to length *trunc*. The underscore method assumes that *flen*, *r* and *trunc*
    are at least 1, and does not support aliasing. Uses binary splitting.

Zeta function
-------------------------------------------------------------------------------

.. function:: void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t s, const arb_t a, int deflate, long n, long prec)

    Sets *res* to the Hurwitz zeta function `\zeta(s,a)` where `s` a power
    series and `a` is a constant, truncated to length *n*.
    To evaluate the usual Riemann zeta function, set `a = 1`.

    If *deflate* is nonzero, evaluates `\zeta(s,a) + 1/(1-s)`, which
    is well-defined as a limit when the constant term of `s` is 1.
    In particular, expanding `\zeta(s,a) + 1/(1-s)` with `s = 1+x`
    gives the Stieltjes constants

    .. math ::

        \sum_{k=0}^{n-1} \frac{(-1)^k}{k!} \gamma_k(a) x^k.

    If `a = 1`, this implementation uses the reflection formula if the midpoint
    of the constant term of `s` is negative.

.. function:: void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    Sets *res* to the series expansion of the Riemann-Siegel theta
    function

    .. math ::

        \theta(h) = \arg \left(\Gamma\left(\frac{2ih+1}{4}\right)\right) - \frac{\log \pi}{2} h

    where the argument of the gamma function is chosen continuously
    as the imaginary part of the log gamma function.

    The underscore method does not support aliasing of the input
    and output arrays, and requires that the lengths are greater
    than zero.

.. function:: void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, long hlen, long n, long prec)

.. function:: void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    Sets *res* to the series expansion of the Riemann-Siegel Z-function

    .. math ::

        Z(h) = e^{i\theta(h)} \zeta(1/2+ih).

    The zeros of the Z-function on the real line precisely
    correspond to the imaginary parts of the zeros of
    the Riemann zeta function on the critical line.

    The underscore method supports aliasing of the input
    and output arrays, and requires that the lengths are greater
    than zero.

Root-finding
-------------------------------------------------------------------------------

.. function:: void _arb_poly_newton_convergence_factor(arf_t convergence_factor, arb_srcptr poly, long len, const arb_t convergence_interval, long prec)

    Given an interval `I` specified by *convergence_interval*, evaluates a bound
    for `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    where `f` is the polynomial defined by the coefficients *{poly, len}*.
    The bound is obtained by evaluating `f'(I)` and `f''(I)` directly.
    If `f` has large coefficients, `I` must be extremely precise in order to
    get a finite factor.

.. function:: int _arb_poly_newton_step(arb_t xnew, arb_srcptr poly, long len, const arb_t x, const arb_t convergence_interval, const arf_t convergence_factor, long prec)

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

.. function:: void _arb_poly_newton_refine_root(arb_t r, arb_srcptr poly, long len, const arb_t start, const arb_t convergence_interval, const arf_t convergence_factor, long eval_extra_prec, long prec)

    Refines a precise estimate of a polynomial root to high precision
    by performing several Newton steps, using nearly optimally
    chosen doubling precision steps.

    The inputs are defined as for *_arb_poly_newton_step*, except for
    the precision parameters: *prec* is the target accuracy and
    *eval_extra_prec* is the estimated number of guard bits that need
    to be added to evaluate the polynomial accurately close to the root
    (typically, if the polynomial has large coefficients of alternating
    signs, this needs to be approximately the bit size of the coefficients).

Other special polynomials
-------------------------------------------------------------------------------

.. function:: void _arb_poly_swinnerton_dyer_ui(arb_ptr poly, ulong n, long trunc, long prec)

.. function:: void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, long prec)

    Computes the Swinnerton-Dyer polynomial `S_n`, which has degree `2^n`
    and is the rational minimal polynomial of the sum
    of the square roots of the first *n* prime numbers.

    If *prec* is set to zero, a precision is chosen automatically such
    that :func:`arb_poly_get_unique_fmpz_poly` should be successful.
    Otherwise a working precision of *prec* bits is used.

    The underscore version accepts an additional *trunc* parameter. Even
    when computing a truncated polynomial, the array *poly* must have room for
    `2^n + 1` coefficients, used as temporary space.

