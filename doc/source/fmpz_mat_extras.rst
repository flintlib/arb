.. _fmpz_mat_extras:

**fmpz_mat_extras.h** -- extra methods for FLINT integer matrices
===============================================================================

This module implements a few utility methods for the FLINT
multiprecision integer matrix type (*fmpz_mat_t*).
It is mainly intended for internal use.

Convenience methods
-------------------------------------------------------------------------------

.. function:: int fmpz_mat_is_hollow(const fmpz_mat_t mat);

    Returns a non-zero value if all entries on the diagonal of *mat* are zero,
    and otherwise returns zero.

.. function:: int fmpz_mat_is_nonnegative(const fmpz_mat_t mat);

    Returns a non-zero value if all entries of *mat* are non-negative,
    and otherwise returns zero.

.. function:: int fmpz_mat_is_lower_triangular(const fmpz_mat_t mat);

    Returns a non-zero value if every entry that is not in the
    lower triangular region of *mat* is zero, and otherwise returns zero.
    The entry at row `i` and column `j` is in the lower triangular region
    if `i >= j`.

.. function:: void fmpz_mat_entrywise_not_is_zero(fmpz_mat_t dest, const fmpz_mat_t src);

    Sets each entry in *dest* to zero or one according to whether
    or not the corresponding entry in *src* is zero.

.. function:: void fmpz_mat_transitive_closure(fmpz_mat_t dest, const fmpz_mat_t src);

    Sets each entry of *dest* to a non-zero value or zero depending
    on whether or not the entry is in the transitive closure of *src*.
    Letting `A` denote the entrywise absolute value of *src*,
    the entry of *dest* at row `i` and column `j` is set to a non-zero
    value if there exists a positive `k` such that the entry in `A^k`
    at row `i` and column `j` is non-zero.
    The *src* matrix must be square.

    Graph-theoretically, this function can be interpreted as computing
    the reachability of `j` from `i` in a graph encoded by *src*.

.. function:: void fmpz_mat_entrywise_nilpotence_degree(fmpz_mat_t B, const fmpz_mat_t A);

    Sets the entry in *B* at row `i` and column `j`
    to the smallest non-negative integer `k` such that the entry in `A^N`
    at row `i` and column `j` is zero for all `N >= k`.
    If no such `k` exists, then the entry in *B* at that row and column
    is set to a negative value.

    If the matrix *B* computed by this function is non-negative
    then this means that *A* is nilpotent.

    This function is implemented only for non-negative square matrices.

    Graph-theoretically, this function can be interpreted as computing
    one plus the longest walk length from `i` to `j` in a graph encoded
    by *src*, with negative values representing infinite walk lengths
    and values of zero or one indicating that no walk exists.
