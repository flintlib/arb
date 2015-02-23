.. _arb-mat:

**arb_mat.h** -- matrices over the real numbers
===============================================================================

An :type:`arb_mat_t` represents a dense matrix over the real numbers,
implemented as an array of entries of type :type:`arb_struct`.

The dimension (number of rows and columns) of a matrix is fixed at
initialization, and the user must ensure that inputs and outputs to
an operation have compatible dimensions. The number of rows or columns
in a matrix can be zero.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: arb_mat_struct

.. type:: arb_mat_t

    Contains a pointer to a flat array of the entries (entries), an array of
    pointers to the start of each row (rows), and the number of rows (r)
    and columns (c).

    An *arb_mat_t* is defined as an array of length one of type
    *arb_mat_struct*, permitting an *arb_mat_t* to
    be passed by reference.

.. macro:: arb_mat_entry(mat, i, j)

    Macro giving a pointer to the entry at row *i* and column *j*.

.. macro:: arb_mat_nrows(mat)

    Returns the number of rows of the matrix.

.. macro:: arb_mat_ncols(mat)

    Returns the number of columns of the matrix.


Memory management
-------------------------------------------------------------------------------

.. function:: void arb_mat_init(arb_mat_t mat, long r, long c)

    Initializes the matrix, setting it to the zero matrix with *r* rows
    and *c* columns.

.. function:: void arb_mat_clear(arb_mat_t mat)

    Clears the matrix, deallocating all entries.


Conversions
-------------------------------------------------------------------------------

.. function:: void arb_mat_set(arb_mat_t dest, const arb_mat_t src)

.. function:: void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src)

.. function:: void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, long prec)

    Sets *dest* to *src*. The operands must have identical dimensions.

Random generation
-------------------------------------------------------------------------------

.. function:: void arb_mat_randtest(arb_mat_t mat, flint_rand_t state, long prec, long mag_bits)

    Sets *mat* to a random matrix with up to *prec* bits of precision
    and with exponents of width up to *mag_bits*.

Input and output
-------------------------------------------------------------------------------

.. function:: void arb_mat_printd(const arb_mat_t mat, long digits)

    Prints each entry in the matrix with the specified number of decimal digits.

Comparisons
-------------------------------------------------------------------------------

.. function:: int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions
    and identical entries.

.. function:: int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions
    and each entry in *mat1* overlaps with the corresponding entry in *mat2*.

.. function:: int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2)

.. function:: int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2)

.. function:: int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions and each entry
    in *mat2* is contained in the corresponding entry in *mat1*.


Special matrices
-------------------------------------------------------------------------------

.. function:: void arb_mat_zero(arb_mat_t mat)

    Sets all entries in mat to zero.

.. function:: void arb_mat_one(arb_mat_t mat)

    Sets the entries on the main diagonal to ones,
    and all other entries to zero.


Norms
-------------------------------------------------------------------------------

.. function:: void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A)

    Sets *b* to an upper bound for the infinity norm (i.e. the largest
    absolute value row sum) of *A*.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void arb_mat_neg(arb_mat_t dest, const arb_mat_t src)

    Sets *dest* to the exact negation of *src*. The operands must have
    the same dimensions.

.. function:: void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)

    Sets res to the sum of *mat1* and *mat2*. The operands must have the same dimensions.

.. function:: void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)

    Sets *res* to the difference of *mat1* and *mat2*. The operands must have
    the same dimensions.

.. function:: void arb_mat_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)

.. function:: void arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)

.. function:: void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)

    Sets *res* to the matrix product of *mat1* and *mat2*. The operands must have
    compatible dimensions for matrix multiplication.

    The *threaded* version splits the computation
    over the number of threads returned by *flint_get_num_threads()*.
    The default version automatically calls the *threaded* version
    if the matrices are sufficiently large and more than one thread
    can be used.

.. function:: void arb_mat_pow_ui(arb_mat_t res, const arb_mat_t mat, ulong exp, long prec)

    Sets *res* to *mat* raised to the power *exp*. Requires that *mat*
    is a square matrix.


Scalar arithmetic
-------------------------------------------------------------------------------

.. function:: void arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, long c)

    Sets *B* to *A* multiplied by `2^c`.

.. function:: void arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)

.. function:: void arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)

.. function:: void arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)

    Sets *B* to `B + A \times c`.

.. function:: void arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)

.. function:: void arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)

.. function:: void arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)

    Sets *B* to `A \times c`.

.. function:: void arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, long c, long prec)

.. function:: void arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)

.. function:: void arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)

    Sets *B* to `A / c`.


Gaussian elimination and solving
-------------------------------------------------------------------------------

.. function:: int arb_mat_lu(long * perm, arb_mat_t LU, const arb_mat_t A, long prec)

    Given an `n \times n` matrix `A`, computes an LU decomposition `PLU = A`
    using Gaussian elimination with partial pivoting.
    The input and output matrices can be the same, performing the
    decomposition in-place.

    Entry `i` in the permutation vector perm is set to the row index in
    the input matrix corresponding to row `i` in the output matrix.

    The algorithm succeeds and returns nonzero if it can find `n` invertible
    (i.e. not containing zero) pivot entries. This guarantees that the matrix
    is invertible.

    The algorithm fails and returns zero, leaving the entries in `P` and `LU`
    undefined, if it cannot find `n` invertible pivot elements.
    In this case, either the matrix is singular, the input matrix was
    computed to insufficient precision, or the LU decomposition was
    attempted at insufficient precision.

.. function:: void arb_mat_solve_lu_precomp(arb_mat_t X, const long * perm, const arb_mat_t LU, const arb_mat_t B, long prec)

    Solves `AX = B` given the precomputed nonsingular LU decomposition `A = PLU`.
    The matrices `X` and `B` are allowed to be aliased with each other,
    but `X` is not allowed to be aliased with `LU`.

.. function:: int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)

    Solves `AX = B` where `A` is a nonsingular `n \times n` matrix
    and `X` and `B` are `n \times m` matrices, using LU decomposition.

    If `m > 0` and `A` cannot be inverted numerically (indicating either that
    `A` is singular or that the precision is insufficient), the values in the
    output matrix are left undefined and zero is returned. A nonzero return
    value guarantees that `A` is invertible and that the exact solution
    matrix is contained in the output.

.. function:: int arb_mat_inv(arb_mat_t X, const arb_mat_t A, long prec)

    Sets `X = A^{-1}` where `A` is a square matrix, computed by solving
    the system `AX = I`.

    If `A` cannot be inverted numerically (indicating either that
    `A` is singular or that the precision is insufficient), the values in the
    output matrix are left undefined and zero is returned.
    A nonzero return value guarantees that the matrix is invertible
    and that the exact inverse is contained in the output.

.. function:: void arb_mat_det(arb_t det, const arb_mat_t A, long prec)

    Computes the determinant of the matrix, using Gaussian elimination
    with partial pivoting. If at some point an invertible pivot element
    cannot be found, the elimination is stopped and the magnitude of the
    determinant of the remaining submatrix is bounded using
    Hadamard's inequality.

Characteristic polynomial
-------------------------------------------------------------------------------

.. function:: void _arb_mat_charpoly(arb_ptr cp, const arb_mat_t mat, long prec)

.. function:: void arb_mat_charpoly(arb_poly_t cp, const arb_mat_t mat, long prec)

    Sets *cp* to the characteristic polynomial of *mat* which must be
    a square matrix. If the matrix has *n* rows, the underscore method
    requires space for `n + 1` output coefficients.
    Employs a division-free algorithm using `O(n^4)` operations.

Special functions
-------------------------------------------------------------------------------

.. function:: void arb_mat_exp(arb_mat_t B, const arb_mat_t A, long prec)

    Sets *B* to the exponential of the matrix *A*, defined by the Taylor series

    .. math ::

        \exp(A) = \sum_{k=0}^{\infty} \frac{A^k}{k!}.

    The function is evaluated as `\exp(A/2^r)^{2^r}`, where `r` is chosen
    to give rapid convergence of the Taylor series. The series is
    evaluated using rectangular splitting.
    If `\|A/2^r\| \le c` and `N \ge 2c`, we bound the entrywise error
    when truncating the Taylor series before term `N` by `2 c^N / N!`.


