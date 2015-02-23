.. _acb-mat:

**acb_mat.h** -- matrices over the complex numbers
===============================================================================

An :type:`acb_mat_t` represents a dense matrix over the complex numbers,
implemented as an array of entries of type :type:`acb_struct`.

The dimension (number of rows and columns) of a matrix is fixed at
initialization, and the user must ensure that inputs and outputs to
an operation have compatible dimensions. The number of rows or columns
in a matrix can be zero.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: acb_mat_struct

.. type:: acb_mat_t

    Contains a pointer to a flat array of the entries (entries), an array of
    pointers to the start of each row (rows), and the number of rows (r)
    and columns (c).

    An *acb_mat_t* is defined as an array of length one of type
    *acb_mat_struct*, permitting an *acb_mat_t* to
    be passed by reference.

.. macro:: acb_mat_entry(mat, i, j)

    Macro giving a pointer to the entry at row *i* and column *j*.

.. macro:: acb_mat_nrows(mat)

    Returns the number of rows of the matrix.

.. macro:: acb_mat_ncols(mat)

    Returns the number of columns of the matrix.


Memory management
-------------------------------------------------------------------------------

.. function:: void acb_mat_init(acb_mat_t mat, long r, long c)

    Initializes the matrix, setting it to the zero matrix with *r* rows
    and *c* columns.

.. function:: void acb_mat_clear(acb_mat_t mat)

    Clears the matrix, deallocating all entries.


Conversions
-------------------------------------------------------------------------------

.. function:: void acb_mat_set(acb_mat_t dest, const acb_mat_t src)

.. function:: void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src)

.. function:: void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, long prec)

    Sets *dest* to *src*. The operands must have identical dimensions.

Random generation
-------------------------------------------------------------------------------

.. function:: void acb_mat_randtest(acb_mat_t mat, flint_rand_t state, long prec, long mag_bits)

    Sets *mat* to a random matrix with up to *prec* bits of precision
    and with exponents of width up to *mag_bits*.

Input and output
-------------------------------------------------------------------------------

.. function:: void acb_mat_printd(const acb_mat_t mat, long digits)

    Prints each entry in the matrix with the specified number of decimal digits.

Comparisons
-------------------------------------------------------------------------------

.. function:: int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions
    and identical entries.

.. function:: int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions
    and each entry in *mat1* overlaps with the corresponding entry in *mat2*.

.. function:: int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2)

.. function:: int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2)

.. function:: int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions and each entry
    in *mat2* is contained in the corresponding entry in *mat1*.

.. function:: int acb_mat_is_real(const acb_mat_t mat)

    Returns nonzero iff all entries in *mat* have zero imaginary part.

Special matrices
-------------------------------------------------------------------------------

.. function:: void acb_mat_zero(acb_mat_t mat)

    Sets all entries in mat to zero.

.. function:: void acb_mat_one(acb_mat_t mat)

    Sets the entries on the main diagonal to ones,
    and all other entries to zero.


Norms
-------------------------------------------------------------------------------

.. function:: void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A)

    Sets *b* to an upper bound for the infinity norm (i.e. the largest
    absolute value row sum) of *A*.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void acb_mat_neg(acb_mat_t dest, const acb_mat_t src)

    Sets *dest* to the exact negation of *src*. The operands must have
    the same dimensions.

.. function:: void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)

    Sets res to the sum of *mat1* and *mat2*. The operands must have the same dimensions.

.. function:: void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)

    Sets *res* to the difference of *mat1* and *mat2*. The operands must have
    the same dimensions.

.. function:: void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)

    Sets *res* to the matrix product of *mat1* and *mat2*. The operands must have
    compatible dimensions for matrix multiplication.

.. function:: void acb_mat_pow_ui(acb_mat_t res, const acb_mat_t mat, ulong exp, long prec)

    Sets *res* to *mat* raised to the power *exp*. Requires that *mat*
    is a square matrix.


Scalar arithmetic
-------------------------------------------------------------------------------

.. function:: void acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, long c)

    Sets *B* to *A* multiplied by `2^c`.

.. function:: void acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)

.. function:: void acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)

.. function:: void acb_mat_scalar_addmul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, long prec)

.. function:: void acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)

    Sets *B* to `B + A \times c`.

.. function:: void acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)

.. function:: void acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)

.. function:: void acb_mat_scalar_mul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, long prec)

.. function:: void acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)

    Sets *B* to `A \times c`.

.. function:: void acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, long c, long prec)

.. function:: void acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)

.. function:: void acb_mat_scalar_div_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, long prec)

.. function:: void acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)

    Sets *B* to `A / c`.


Gaussian elimination and solving
-------------------------------------------------------------------------------

.. function:: int acb_mat_lu(long * perm, acb_mat_t LU, const acb_mat_t A, long prec)

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

.. function:: void acb_mat_solve_lu_precomp(acb_mat_t X, const long * perm, const acb_mat_t LU, const acb_mat_t B, long prec)

    Solves `AX = B` given the precomputed nonsingular LU decomposition `A = PLU`.
    The matrices `X` and `B` are allowed to be aliased with each other,
    but `X` is not allowed to be aliased with `LU`.

.. function:: int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)

    Solves `AX = B` where `A` is a nonsingular `n \times n` matrix
    and `X` and `B` are `n \times m` matrices, using LU decomposition.

    If `m > 0` and `A` cannot be inverted numerically (indicating either that
    `A` is singular or that the precision is insufficient), the values in the
    output matrix are left undefined and zero is returned. A nonzero return
    value guarantees that `A` is invertible and that the exact solution
    matrix is contained in the output.

.. function:: int acb_mat_inv(acb_mat_t X, const acb_mat_t A, long prec)

    Sets `X = A^{-1}` where `A` is a square matrix, computed by solving
    the system `AX = I`.

    If `A` cannot be inverted numerically (indicating either that
    `A` is singular or that the precision is insufficient), the values in the
    output matrix are left undefined and zero is returned.
    A nonzero return value guarantees that the matrix is invertible
    and that the exact inverse is contained in the output.

.. function:: void acb_mat_det(acb_t det, const acb_mat_t A, long prec)

    Computes the determinant of the matrix, using Gaussian elimination
    with partial pivoting. If at some point an invertible pivot element
    cannot be found, the elimination is stopped and the magnitude of the
    determinant of the remaining submatrix is bounded using
    Hadamard's inequality.

Characteristic polynomial
-------------------------------------------------------------------------------

.. function:: void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, long prec)

.. function:: void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, long prec)

    Sets *cp* to the characteristic polynomial of *mat* which must be
    a square matrix. If the matrix has *n* rows, the underscore method
    requires space for `n + 1` output coefficients.
    Employs a division-free algorithm using `O(n^4)` operations.

Special functions
-------------------------------------------------------------------------------

.. function:: void acb_mat_exp(acb_mat_t B, const acb_mat_t A, long prec)

    Sets *B* to the exponential of the matrix *A*, defined by the Taylor series

    .. math ::

        \exp(A) = \sum_{k=0}^{\infty} \frac{A^k}{k!}.

    The function is evaluated as `\exp(A/2^r)^{2^r}`, where `r` is chosen
    to give rapid convergence of the Taylor series. The series is
    evaluated using rectangular splitting.
    If `\|A/2^r\| \le c` and `N \ge 2c`, we bound the entrywise error
    when truncating the Taylor series before term `N` by `2 c^N / N!`.

