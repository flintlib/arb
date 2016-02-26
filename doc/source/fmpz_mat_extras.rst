.. _fmpz_mat_extras:

**fmpz_mat_extras.h** -- extra methods for FLINT integer matrices
===============================================================================

This module implements a few utility methods for the FLINT
multiprecision integer matrix type (*fmpz_mat_t*).
It is mainly intended for internal use.


Comparison
-------------------------------------------------------------------------------

.. function:: int fmpz_mat_is_hollow(const fmpz_mat_t mat)

    Returns a non-zero value if all entries on the diagonal of *mat* are zero,
    and otherwise returns zero.

.. function:: int fmpz_mat_is_nonnegative(const fmpz_mat_t mat)

    Returns a non-zero value if all entries of *mat* are non-negative,
    and otherwise returns zero.

.. function:: int fmpz_mat_is_lower_triangular(const fmpz_mat_t mat)

    Returns a non-zero value if every entry that is not in the
    lower triangular region of *mat* is zero, and otherwise returns zero.
    The entry at row `i` and column `j` is in the lower triangular region
    if `i \ge j`.

.. function:: void fmpz_mat_entrywise_not_is_zero(fmpz_mat_t dest, const fmpz_mat_t src)

    Sets each entry in *dest* to zero or one according to whether
    or not the corresponding entry in *src* is zero.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpz_mat_add_ui_entrywise(fmpz_mat_t B, const fmpz_mat_t A, ulong x)

    Sets `B_{ij}` to `A_{ij} + x` where `x` is an ulong.

.. function:: void fmpz_mat_sub_ui_entrywise(fmpz_mat_t B, const fmpz_mat_t A, ulong x)

    Sets `B_{ij}` to `A_{ij} - x` where `x` is an ulong.


Graph theory
-------------------------------------------------------------------------------

.. function:: void fmpz_mat_transitive_closure(fmpz_mat_t B, const fmpz_mat_t A)

    Given the adjacency matrix *A* of an unweighted directed graph `G`,
    *B* is set to the adjacency matrix of the transitive closure of `G`.
    In particular, `B_{st}` is set to a non-zero value if there is a non-empty
    walk from `s` to `t` in `G`, and otherwise `B_{st}` is set to zero.

.. function:: void fmpz_mat_unweighted_all_pairs_longest_walk(fmpz_mat_t B, const fmpz_mat_t A)

    Given the adjacency matrix *A* of an unweighted directed graph `G`,
    `B_{st}` is set to the length of the longest walk from `s` to `t` in `G`.
    Empty walks are considered, so for each vertex `s` there is a zero-length
    walk from `s` to `s`.  If `t` is unreachable from `s` then `B_{st}`
    is set to the special value `-1`.  If `G` contains a cycle then arbitrarily
    long walks may exist; if the length of the longest walk from `s` to `t`
    is unbounded then `B_{st}` is set to the special value `-2`.

    This function can help quantify entrywise errors in a truncated evaluation
    of a matrix power series.  If *A* is an indictor matrix with the same
    sparsity pattern as a matrix `M`, and if `B_{ij}` does not take
    the special value `-2`, then the tail
    `\left[ \sum_{k=N}^\infty a_k M^k \right]_{ij}`
    vanishes when `N > B_{ij}`.
