fmpcb_poly.h -- polynomials over the complex numbers
===============================================================================

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

    Returns the length of the polynomial.

.. function:: void fmpcb_poly_zero(fmpcb_poly_t poly)

    Sets *poly* to the zero polynomial.

.. function:: void fmpcb_poly_one(fmpcb_poly_t poly)

    Sets *poly* to the constant polynomial 1.

.. function:: void fmpcb_poly_set(fmpcb_poly_t dest, const fmpcb_poly_t src)

    Sets *dest* to a copy of *src*.

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

.. function:: int fmpcb_poly_contains_fmpq_poly(const fmpcb_poly_t poly1, const fmpq_poly_t poly2)

.. function:: int fmpcb_poly_contains(const fmpcb_poly_t poly1, const fmpcb_poly_t poly2)

    Returns nonzero iff *poly2* is contained in *poly1*.

.. function:: int _fmpcb_poly_overlaps(const fmpcb_struct * poly1, long len1,
        const fmpcb_struct * poly2, long len2)

.. function:: int fmpcb_poly_overlaps(const fmpcb_poly_t poly1, const fmpcb_poly_t poly2)

    Returns nonzero iff *poly1* overlaps with *poly2*. The underscore
    function requires that *len1* ist at least as large as *len2*.


Conversions
-------------------------------------------------------------------------------

.. function:: void fmpcb_poly_set_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re)

.. function:: void fmpcb_poly_set2_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re, const fmprb_poly_t im)

.. function:: void fmpcb_poly_set_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, long prec)

.. function:: void fmpcb_poly_set2_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec)

    Sets *poly* to the given real polynomial *re* plus the polynomial *im*
    multiplied by the imaginary unit.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_add(fmpcb_struct * res, const fmpcb_struct * poly1, long len1, const fmpcb_struct * poly2, long len2, long prec)

.. function:: void fmpcb_poly_add(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long prec)

.. function:: void _fmpcb_poly_mullow_classical(fmpcb_struct * res, const fmpcb_struct * poly1, long len1, const fmpcb_struct * poly2, long len2, long n, long prec)

.. function:: void fmpcb_poly_mullow_classical(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long n, long prec)

.. function:: void _fmpcb_poly_mullow_transpose(fmpcb_struct * res, const fmpcb_struct * poly1, long len1, const fmpcb_struct * poly2, long len2, long n, long prec)

.. function:: void fmpcb_poly_mullow_transpose(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long n, long prec)

.. function:: void _fmpcb_poly_mullow(fmpcb_struct * res, const fmpcb_struct * poly1, long len1, const fmpcb_struct * poly2, long len2, long n, long prec)

.. function:: void fmpcb_poly_mullow(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long n, long prec)

.. function:: void _fmpcb_poly_mul(fmpcb_struct * C, const fmpcb_struct * A, long lenA, const fmpcb_struct * B, long lenB, long prec)

.. function:: void fmpcb_poly_mul(fmpcb_poly_t res, const fmpcb_poly_t poly1, const fmpcb_poly_t poly2, long prec)

.. function:: void _fmpcb_poly_inv_series(fmpcb_struct * Qinv, const fmpcb_struct * Q, long len, long prec)

.. function:: void fmpcb_poly_inv_series(fmpcb_poly_t Qinv, const fmpcb_poly_t Q, long n, long prec)


Evaluation
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_evaluate(fmpcb_t res, const fmpcb_struct * f, long len, const fmpcb_t a, long prec)

.. function:: void fmpcb_poly_evaluate(fmpcb_t res, const fmpcb_poly_t f, const fmpcb_t a, long prec)

    Evaluates the polynomial using Horner's rule.


Derivatives
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_derivative(fmpcb_struct * res, const fmpcb_struct * poly, long len, long prec)

.. function:: void fmpcb_poly_derivative(fmpcb_poly_t res, const fmpcb_poly_t poly, long prec)

    Sets *res* to the derivative of *poly*.


Root-finding
-------------------------------------------------------------------------------

.. function:: void _fmpcb_poly_root_inclusion(fmpcb_t r, const fmpcb_t m, const fmpcb_struct * poly, const fmpcb_struct * polyder, long len, long prec);

    Given any complex number `m`, and a nonconstant polynomial `f` and its
    derivative `f'`, sets *r* to an interval centered on `m` that is
    guaranteed to contain at least one root of `f`.
    Such an interval is obtained by taking a ball of radius `|f(m)/f'(m)| n`
    where `n` is the degree of `f`.

.. function:: long _fmpcb_poly_validate_roots(fmpcb_struct * roots, const fmpcb_struct * poly, long len, long prec)

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

.. function:: void _fmpcb_poly_refine_roots_durand_kerner(fmpcb_struct * roots, const fmpcb_struct * poly, long len, long prec)

    Refines the given roots simultaneously using a single iteration
    of the Durand-Kerner method. The radius of each root is set to an
    approximation of the correction, giving a rough estimate of its error (not
    a rigorous bound).

.. function:: long _fmpcb_poly_find_roots(fmpcb_struct * roots, const fmpcb_struct * poly, const fmpcb_struct * initial, long len, long maxiter, long prec)

.. function:: long fmpcb_poly_find_roots(fmpcb_struct * roots, const fmpcb_poly_t poly, const fmpcb_struct * initial, long maxiter, long prec)

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

