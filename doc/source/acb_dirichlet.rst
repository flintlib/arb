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

Multiplicative group modulo *q*
-------------------------------------------------------------------------------

Working with Dirichlet characters mod *q* consists mainly
in going from residue classes mod *q* to exponents on a set
of generators of the group.

This implementation relies on the Conrey numbering scheme
introduced in the LMFDB.
We call *number* a residue class modulo *q*, and *index* the
corresponding vector of exponents of Conrey generators.

Going from an *index* to the corresponding *number* is a cheap
operation while the converse requires computing discrete
logarithms.

.. type:: acb_dirichlet_group_struct

.. type:: acb_dirichlet_group_t

    Represents the group of Dirichlet characters mod *q*.

    An *acb_dirichlet_group_t* is defined as an array of *acb_dirichlet_group_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q)

    Initializes *G* to the group of Dirichlet characters mod *q*.

    This method computes a canonical decomposition of *G* in terms of cyclic
    groups, which are the mod `p^e` subgroups for `p^e\|q`.
    In particular *G* contains:

    - the number *num* of components

    - the generators

    - the exponent *expo* of the group

    It does *not* automatically precompute lookup tables
    of discrete logarithms or numerical roots of unity, and can therefore
    safely be called even with large *q*.

    For implementation reasons, the largest prime factor of *q* must not
    exceed `10^{12}` (an abort will be raised). This restriction could
    be removed in the future.

.. function:: void acb_dirichlet_group_clear(acb_dirichlet_group_t G)

    Clears *G*.

.. function:: void acb_dirichlet_group_dlog_precompute(acb_dirichlet_group_t G, ulong num)

    Precompute decomposition and tables for discrete log computations in *G*,
    so as to minimize the complexity of *num* calls to discrete logarithms.

    If *num* gets very large, the entire group may be indexed.

Conrey index
...............................................................................

.. type:: acb_dirichlet_conrey_struct

.. type:: acb_dirichlet_conrey_t

    Represents elements of the unit group mod *q*, keeping both the
    *number* (residue class) and *index* (exponents on the group
    generators).

.. function:: void acb_dirichlet_conrey_log(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G, ulong m)

    Sets *x* to the element of number *m*, computing its index using discrete
    logarithm in *G*.

.. function:: ulong acb_dirichlet_conrey_exp(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G)

    Compute the reverse operation.

.. function:: void acb_dirichlet_conrey_one(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G)

    Sets *x* to the *number* `1\in G`, having *index* `[0,\dots 0]`.

.. function:: int acb_dirichlet_conrey_next(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G)

    This function allows to iterate on the elements of *G* looping on the *index*.
    It produces elements in seemingly random *number* order. The return value
    is the index of the last updated exponent of *x*, or *G->num* if the last
    element has been reached.

Dirichlet characters
-------------------------------------------------------------------------------

Dirichlet characters take value in a finite cyclic group of roots of unity plus zero.

When evaluation functions return a *ulong*, this number corresponds to the
power of a primitive root of unity, the special value *ACB_DIRICHLET_CHI_NULL*
encoding the zero value.

The Conrey numbering scheme makes explicit the mathematical fact that
the group *G* is isomorphic to its dual.

.. function:: ulong acb_dirichlet_ui_pairing(const acb_dirichlet_group_t G, ulong m, ulong n)

.. function:: ulong acb_dirichlet_ui_pairing_conrey(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, const acb_dirichlet_conrey_t b)

    Compute the value of the Dirichlet pairing on numbers *m* and *n*, as
    exponent modulo *G->expo*.
    The second form takes the index *a* and *b*, and does not take discrete
    logarithms.

    The returned value is the numerator of the actual value exponent mod the group exponent *G->expo*.

Character properties
...............................................................................

As a consequence of the Conrey numbering, properties of
characters such that the order, the parity or the conductor are available
at the level of their number, whithout any discrete log computation,
or at the Conrey index level.

.. function:: ulong acb_dirichlet_ui_order(const acb_dirichlet_group_t G, ulong a)

.. function:: ulong acb_dirichlet_conrey_order(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

   Compute the order of a mod q.

.. function:: ulong acb_dirichlet_ui_conductor(const acb_dirichlet_group_t G, ulong a)

.. function:: ulong acb_dirichlet_conrey_conductor(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

   Compute the conductor of a mod q, that is the smallest r dividing q such
   that a corresponds to an element defined modulo r.

.. function:: ulong acb_dirichlet_ui_parity(const acb_dirichlet_group_t G, ulong a)

.. function:: int acb_dirichlet_conrey_parity(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

   Compute the parity of a mod q, which is the parity of the corresponding
   Dirichlet character.

Character type
-------------------------------------------------------------------------------

.. type:: acb_dirichlet_char_struct

.. type:: acb_dirichlet_char_t

    Represents a Dirichlet character. This structure contains various
    useful invariants such as the order, the parity and the conductor of the character.

    An *acb_dirichlet_char_t* is defined as an array of *acb_dirichlet_char_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dirichlet_char_init(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);

.. function:: void acb_dirichlet_char_clear(acb_dirichlet_char_t chi);

    Initializes and clear *chi*.

.. function:: void acb_dirichlet_char(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, ulong n);

    Sets *chi* to the Dirichlet character of number *n*, using Conrey numbering scheme.
    This function performs a discrete logarithm in *G*.

.. function:: void acb_dirichlet_char_conrey(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);

    Sets *chi* to the Dirichlet character of Conrey index *x*.

Character properties
...............................................................................

.. function:: ulong acb_dirichlet_char_order(const acb_dirichlet_char_t chi)

.. function:: int acb_dirichlet_char_parity(const acb_dirichlet_char_t chi)

.. function:: ulong acb_dirichlet_char_conductor(const acb_dirichlet_char_t chi)

Character evaluation
-------------------------------------------------------------------------------

The image of a Dirichlet character is a finite cyclic group. Dirichlet
character evaluations are either exponents in this group, or an *acb_t* root of
unity.

.. function:: void acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n, slong prec)

    Sets *res* to `\chi(n)`, the value of the Dirichlet character *chi*
    at the integer *n*.

    There are no restrictions on *n*.

Gauss and Jacobi sums
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_gauss_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_jacobi_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1,  const acb_dirichlet_char_t chi2, slong prec)

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


Implementation notes
-------------------------------------------------------------------------------

The current implementation introduces a *char* type which contains a *conrey*
index plus additional information which

- make evaluation of a single character a bit faster

- have some initialization cost.

Even if it is straiforward to convert a *conrey* index to the
corresponding *char*, looping is faster at the
level of conrey representation. Things can be improved on this aspect
but it makes code more intricate.
