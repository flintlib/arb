.. _acb-dirichlet:

**acb_dirichlet.h** -- Dirichlet L-functions, zeta functions, and related functions
===================================================================================

Warning: the interfaces in this module are experimental and may change
without notice.

This module will eventually allow working with Dirichlet L-functions and
possibly slightly more general Dirichlet series. At the moment, it contains
nothing interesting.

The code in other modules for computing the Riemann zeta function,
Hurwitz zeta function and polylogarithm will possibly be migrated to this
module in the future.

A Dirichlet L-function is the analytic continuation of an L-series

.. math ::

    L(s,\chi) = \sum_{k=1}^\infty \frac{\chi(k)}{k^s}

where `\chi(k)` is a Dirichlet character.

Dirichlet characters
-------------------------------------------------------------------------------

.. type:: acb_dirichlet_struct

.. type:: acb_dirichlet_t

    Represents the group of Dirichlet characters mod *q*.

    An *acb_dirichlet_t* is defined as an array of *acb_dirichlet_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q)

    Initializes *G* to the group of Dirichlet characters mod *q*.

    This method computes the prime factorization of *q* and other useful
    invariants. It does *not* automatically precompute lookup tables
    of discrete logarithms or numerical roots of unity, and can therefore
    safely be called even with large *q*.
    For implementation reasons, the largest prime factor of *q* must not
    exceed `10^{12}` (an abort will be raised). This restriction could
    be removed in the future.

.. function:: void acb_dirichlet_group_clear(acb_dirichlet_group_t G)

    Clears *G*.

.. function:: void acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec)

    Sets *res* to `\chi_m(n)`, the value of the Dirichlet character
    of index *m* evaluated at the integer *n*.

    Requires that *m* is a valid index, that is, `1 \le m \le q` and *m* is
    coprime to *q*. There are no restrictions on *n*.

Euler products
-------------------------------------------------------------------------------

.. function:: void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s, const signed char * chi, int mod, int reciprocal, slong prec)

    Sets *res* to `L(s,\chi)` where `\chi` is a real Dirichlet character
    given by the explicit list *chi* of character values at
    0, 1, ..., *mod* - 1. If *reciprocal* is set, computes `1 / L(s,\chi)`
    (this is faster if the reciprocal can be used directly).

    This function uses the Euler product, and is only intended for use when
    *s* is large. An error bound is computed via :func:`mag_hurwitz_zeta_uiui`.
    Since

    .. math ::

        \frac{1}{L(s,\chi)} = \prod_{p} \left(1 - \frac{\chi(p)}{p^s}\right)
                = \sum_{k=1}^{\infty} \frac{\mu(k)\chi(k)}{k^s}

    and the truncated product gives all smooth-index terms in the series, we have

    .. math ::

        \left|\prod_{p < N} \left(1 - \frac{\chi(p)}{p^s}\right) - \frac{1}{L(s,\chi)}\right|
        \le \sum_{k=N}^{\infty} \frac{1}{k^s} = \zeta(s,N).

Simple functions
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec)

    Sets *res* to the Dirichlet eta function
    `\eta(s) = \sum_{k=1}^{\infty} (-1)^k / k^s = (1-2^{1-s}) \zeta(s)`,
    also known as the alternating zeta function.
    Note that the alternating character `\{1,-1\}` is not itself
    a Dirichlet character.

