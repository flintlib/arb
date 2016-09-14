.. _acb-dirichlet:

**acb_dirichlet.h** -- Dirichlet L-functions, zeta functions, and related functions
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

This module allows working with Dirichlet characters, Dirichlet L-functions,
and related functions.
A Dirichlet L-function is the analytic continuation of an L-series

.. math ::

    L(s,\chi) = \sum_{k=1}^\infty \frac{\chi(k)}{k^s}

where `\chi(k)` is a Dirichlet character.

The code in other modules for computing the Riemann zeta function,
Hurwitz zeta function and polylogarithm will possibly be migrated to this
module in the future.

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

.. function:: void acb_dirichlet_subgroup_init(acb_dirichlet_group_t H, const acb_dirichlet_group_t G, ulong h)

   Given an already computed group *G* mod `q`, initialize its subgroup *H*
   defined mod `h\mid q`. Precomputed discrete logs tables are kept.

.. function:: void acb_dirichlet_group_clear(acb_dirichlet_group_t G)

    Clears *G*.

.. function:: void acb_dirichlet_group_dlog_precompute(acb_dirichlet_group_t G, ulong num)

    Precompute decomposition and tables for discrete log computations in *G*,
    so as to minimize the complexity of *num* calls to discrete logarithms.

    If *num* gets very large, the entire group may be indexed.

Conrey elements
-------------------------------------------------------------------------------

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

    Sets *x* to the next conrey index in *G* with lexicographic ordering.
    The return value
    is the index of the last updated exponent of *x*, or *-1* if the last
    element has been reached.

    This function allows to iterate on the elements of *G* looping on their *index*.
    Note that it produces elements in seemingly random *number* order.

    The following template can be used to loop over all elements *x* in *G*::

        acb_conrey_one(x, G);
        do {
            /* use Conrey element x */
        } while (acb_dirichlet_conrey_next(x, G) >= 0);


.. function:: int acb_dirichlet_conrey_eq(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x, const acb_dirichlet_conrey_t y)

   Return 1 if *x* equals *y*. This function checks both *number* and *index*,
   writing ``(x->n == y->n)`` gives a faster check.

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
   The second form takes the Conrey index *a* and *b*, and does not take discrete
   logarithms.

   The returned value is the numerator of the actual value exponent mod the group exponent *G->expo*.

Character type
-------------------------------------------------------------------------------

.. type:: acb_dirichlet_char_struct

.. type:: acb_dirichlet_char_t

    Represents a Dirichlet character. This structure contains various
    useful invariants such as the order, the parity and the conductor of the character.

    An *acb_dirichlet_char_t* is defined as an array of *acb_dirichlet_char_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dirichlet_char_init(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)

.. function:: void acb_dirichlet_char_clear(acb_dirichlet_char_t chi)

    Initializes and clear *chi*.

.. function:: void acb_dirichlet_char(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, ulong n)

    Sets *chi* to the Dirichlet character of number *n*, using Conrey numbering scheme.
    This function performs a discrete logarithm in *G*.

.. function:: void acb_dirichlet_char_conrey(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

    Sets *chi* to the Dirichlet character of Conrey index *x*.

.. function:: int acb_dirichlet_char_eq(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2)

   Return 1 if *chi1* equals *chi2*.

.. function:: acb_dirichlet_char_is_principal(const acb_dirichlet_char_t chi)

   Return 1 if *chi* is the principal character mod *q*.

Character properties
-------------------------------------------------------------------------------

As a consequence of the Conrey numbering, all these numbers are available at the
level of *number* and *index*, and for *char*.
No discrete log computation is performed.

.. function:: ulong acb_dirichlet_number_primitive(const acb_dirichlet_group_t G)

   Return the number of primitive elements in *G*.

.. function:: ulong acb_dirichlet_ui_conductor(const acb_dirichlet_group_t G, ulong a)

.. function:: ulong acb_dirichlet_conrey_conductor(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

.. function:: ulong acb_dirichlet_char_conductor(const acb_dirichlet_char_t chi)

   Return the *conductor* of `\chi_q(a,\cdot)`, that is the smallest `r` dividing `q`
   such `\chi_q(a,\cdot)` can be obtained as a character mod `r`.
   This number is precomputed for the *char* type.

.. function:: int acb_dirichlet_ui_parity(const acb_dirichlet_group_t G, ulong a)

.. function:: int acb_dirichlet_conrey_parity(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

.. function:: int acb_dirichlet_char_parity(const acb_dirichlet_char_t chi)

   Return the *parity* `\lambda` in `\{0,1\}` of `\chi_q(a,\cdot)`, such that
   `\chi_q(a,-1)=(-1)^\lambda`.
   This number is precomputed for the *char* type.

.. function:: ulong acb_dirichlet_ui_order(const acb_dirichlet_group_t G, ulong a)

.. function:: int acb_dirichlet_conrey_order(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)

.. function:: ulong acb_dirichlet_char_order(const acb_dirichlet_char_t chi)

   Return the order of `\chi_q(a,\cdot)` which is the order of `a\bmod q`.
   This number is precomputed for the *char* type.

.. function:: acb_dirichlet_char_is_real(const acb_dirichlet_char_t chi)

   Return 1 if *chi* is a real character (iff it has order `\leq 2`).

Character evaluation
-------------------------------------------------------------------------------

The image of a Dirichlet character is a finite cyclic group. Dirichlet
character evaluations are either exponents in this group, or an *acb_t* root of
unity.

.. function:: ulong acb_dirichlet_ui_chi_conrey(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const acb_dirichlet_conrey_t x)

.. function:: ulong acb_dirichlet_ui_chi(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n)

   Compute that value `\chi(a)` as the exponent mod the order of `\chi`.

.. function:: void acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n, slong prec)

    Sets *res* to `\chi(n)`, the value of the Dirichlet character *chi*
    at the integer *n*.

    There are no restrictions on *n*.

Roots of unity
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_nth_root(acb_t res, ulong order, slong prec)

   Sets *res* to `\exp(\frac{2i\pi}{\mathrm{order}})` to precision *prec*.

.. function:: void acb_dirichlet_vec_nth_roots(acb_ptr z, slong order, slong prec)

   Compute the vector ``1,z,z^2,\dots z^{\mathrm{order}-1}`` where `z=\exp(\frac{2i\pi}{\mathrm{order}})` to precision *prec*.

   In order to avoid precision loss, this function does not simply compute powers of a primitive root.

.. type:: acb_dirichlet_powers_struct

.. type:: acb_dirichlet_powers_t

   This structure allows to compute *n* powers of a fixed root of unity of order *m*
   using precomputations. Extremal cases are

   - all powers are stored: `O(m)` initialization + storage, `O(n)` evaluations

   - nothing stored: `O(1)` initialization + storage, `O(\log(m)n)` evaluations

   - `k` step decomposition: `O(k m^{\frac1k})` init + storage, `O((k-1)n)` evaluations.

   Currently, only baby-step giant-step decomposition (i.e. `k=2`)
   is implemented, allowing to obtain each power using one multiplication.

.. function:: void acb_dirichlet_powers_init(acb_dirichlet_powers_t t, ulong order, slong num, slong prec)

   Initialize the powers structure for *num* evaluations of powers of the root of unity
   of order *order*.

.. function:: void acb_dirichlet_powers_clear(acb_dirichlet_powers_t t)

   Clears *t*.

.. function:: void acb_dirichlet_power(acb_t z, const acb_dirichlet_powers_t t, ulong n, slong prec)

   Sets *z* to `x^n` where *t* contains precomputed powers of `x`.

Vector evaluation
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_ui_chi_vec(ulong * v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)

   Compute the list of exponent values *v[k]* for `0\leq k < nv`.

.. function:: void acb_dirichlet_chi_vec(acb_ptr v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv, slong prec)

   Compute the *nv* first Dirichlet values.

Character operations
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_conrey_mul(acb_dirichlet_conrey_t c, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, const acb_dirichlet_conrey_t b)

.. function:: void acb_dirichlet_char_mul(acb_dirichlet_char_t chi12, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2)

   Multiply two characters in the same group.

.. function:: void acb_dirichlet_conrey_pow(acb_dirichlet_conrey_t c, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, ulong n)

   Take the power of some character.

Gauss and Jacobi sums
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_gauss_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)

   Compute the Gauss sum

   .. math::

      G_q(a) = \sum_{x \bmod q} \chi_q(a, x)e^{\frac{2i\pi x}q}

.. function:: void acb_dirichlet_jacobi_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1,  const acb_dirichlet_char_t chi2, slong prec)

   Compute the Jacobi sum

   .. math::

      J_q(a,b) = \sum_{x \bmod q} \chi_q(a, x)\chi_q(b, 1-x)

Theta sums
-------------------------------------------------------------------------------

We call *theta series* of a Dirichlet character the quadratic series

.. math::

   \Theta_q(a) = \sum_{n\geq 0} \chi_q(a, n) n^p x^{n^2}

where `p` is the parity of the character `\chi_q(a,\cdot)`.

For `\Re(t)>0` we write `x(t)=\exp(-\frac{\pi}{N}t^2)` and define

.. math::

   \Theta_q(a,t) = \sum_{n\geq 0} \chi_q(a, n) x(t)^{n^2}.

.. function:: void acb_dirichlet_chi_theta_arb(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t t, slong prec)

.. function:: void acb_dirichlet_ui_theta_arb(acb_t res, const acb_dirichlet_group_t G, ulong a, const arb_t t, slong prec)

   Compute the theta series `\Theta_q(a,t)` for real argument `t>0`.
   Beware that if `t<1` the functional equation

   .. math::

      t \theta(a,t) = \epsilon(\chi) \theta\left(\frac1a, \frac1t\right)

   should be used, which is not done automatically (to avoid recomputing the
   Gauss sum).

.. function:: ulong acb_dirichlet_theta_length(ulong q, const arb_t t, slong prec)

   Compute the number of terms to be summed in the theta series of argument *t*
   so that the tail is less than `2^{-\mathrm{prec}}`.

.. function:: void acb_dirichlet_qseries_powers_naive(acb_t res, const arb_t x, int p, const ulong * a, const acb_dirichlet_powers_t z, slong len, slong prec)

.. function:: void acb_dirichlet_qseries_powers_smallorder(acb_t res, const arb_t x, int p, const ulong * a, const acb_dirichlet_powers_t z, slong len, slong prec)

   Compute the series `\sum n^p z^{a_n} x^{n^2}` for exponent list *a*,
   precomputed powers *z* and parity *p* (being 0 or 1).

   The *naive* version sums the series as defined, while the *smallorder*
   variant evaluates the series on the quotient ring by a cyclotomic polynomial
   before evaluating at the root of unity, ignoring its argument *z*.

Discrete Fourier Transforms (DFT)
-------------------------------------------------------------------------------

Let *G* be a finite abelian group, and `\chi` a character of *G*.
For any map `f:G\to\mathbb C`, the discrete fourier transform
`\hat f:\hat G\to \mathbb C` is defined by

.. math::

   \hat f(\chi) = \sum_{x\in G}\chi(x)f(x)

Fast Fourier Transform techniques allow to compute efficiently
all values `\hat f(\chi)`.

For a Dirichlet group `G` modulo `q`, we take advantage
of the Conrey isomorphism `G \to \hat G` to consider the
the Fourier transform on Conrey labels as

.. math::

   g(a) = \sum_{b\bmod q}\chi_q(a,b)f(b)


.. function:: void acb_dirichlet_dft_conrey(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec)

   Compute the DFT of *v* using Conrey indices.
   This function assumes *v* and *w* are vectors
   of size *G->phi_q*, whose values correspond to a lexicographic ordering
   of Conrey indices (as obtained using :func:`acb_dirichlet_conrey_next`).

   For example, if `q=15`, the Conrey elements are stored in following
   order

   ============  =====================
   index [e,f]     number = 7^e 11^f
   ============  =====================
   [0, 0]        1
   [0, 1]        7
   [0, 2]        4
   [0, 3]        13
   [0, 4]        1
   [1, 0]        11
   [1, 1]        2
   [1, 2]        14
   [1, 3]        8
   [1, 4]        11
   ============  =====================

.. function:: void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec)

   Compute the DFT of *v* using Conrey numbers.
   This function assumes *v* and *w* are vectors of size *G->q*.
   All values at index not coprime to *G->q* are ignored.

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

L-functions
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)

    Compute `L(s,\chi)` using decomposition in terms of the Hurwitz zeta function

    .. math::

        L(s,\chi) = q^{-s}\sum_{k=1}^{q-1} \chi(k) \,\zeta\!\left(s,\frac kq\right).

    If `s = 1` and `\chi` is non-principal, the deflated Hurwitz zeta function
    is used to avoid poles.

    This formula is slow for large *q*.

.. function:: void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const acb_dirichlet_group_t G, slong prec)

    Compute all values `L(s,\chi)` for `\chi` mod `q`, by Hurwitz formula and
    discrete Fourier transform.
    *res* is assumed to have length *G->phi_q* and values are stored by lexicographically ordered Conrey
    index. See :func:`acb_dirichlet_dft_conrey`.

Implementation notes
-------------------------------------------------------------------------------

The current implementation introduces a *char* type which contains a *conrey*
index plus additional information which

- makes evaluation of a single character a bit faster

- has some initialization cost.

Even if it is straightforward to convert a *conrey* index to the
corresponding *char*, looping is faster at the
level of Conrey representation. Things can be improved on this aspect
but it makes code more intricate.
