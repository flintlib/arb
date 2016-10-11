.. _acb-dirichlet:

**acb_dirichlet.h** -- Dirichlet L-functions, zeta functions, and related functions
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

This module allows working with values of Dirichlet characters, Dirichlet L-functions,
and related functions. Working with Dirichlet characters is documented in
:ref:`dirichlet`.

A Dirichlet L-function is the analytic continuation of an L-series

.. math ::

    L(s,\chi) = \sum_{k=1}^\infty \frac{\chi(k)}{k^s}

where `\chi(k)` is a Dirichlet character.

The code in other modules for computing the Riemann zeta function,
Hurwitz zeta function and polylogarithm will possibly be migrated to this
module in the future.

Character evaluation
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, slong prec)

    Sets *res* to `\chi(n)`, the value of the Dirichlet character *chi*
    at the integer *n*.

.. function:: void acb_dirichlet_chi_vec(acb_ptr v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv, slong prec)

    Compute the *nv* first Dirichlet values.

.. function:: void acb_dirichlet_pairing(acb_t res, const dirichlet_group_t G, ulong m, ulong n, slong prec)

.. function:: void acb_dirichlet_pairing_char(acb_t res, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b, slong prec)

    Sets *res* to the value of the Dirichlet pairing `\chi(m,n)` at numbers `m` and `n`.
    The second form takes two characters as input.

Roots of unity
-------------------------------------------------------------------------------

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

Gauss and Jacobi sums
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_gauss_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_gauss_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_gauss_sum_order2(acb_t res, const dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_gauss_sum_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_gauss_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_gauss_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, slong prec)

   Sets *res* to the Gauss sum

   .. math::

      G_q(a) = \sum_{x \bmod q} \chi_q(a, x) e^{\frac{2i\pi x}q}

   - the *naive* version computes the sum as defined.

   - the *factor* version writes it as a product of local Gauss sums by chinese
     remainder theorem.

   - the *order2* version assumes *chi* is real and primitive and returns
     `i^p\sqrt q` where `p` is the parity of `\chi`.

   - the *theta* version assumes that *chi* is primitive to obtain the Gauss
     sum by functional equation of the theta series at `t=1`. An abort will be
     raised if the theta series vanishes at `t=1`. Only 4 exceptional
     characters of conductor 300 and 600 are known to have this particularity,
     and none with primepower modulus.

   - the default version automatically combines the above methods.

   - the *ui* version only takes the Conrey number *a* as parameter.

.. function:: void acb_dirichlet_jacobi_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)

.. function:: void acb_dirichlet_jacobi_sum_factor(acb_t res,  const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)

.. function:: void acb_dirichlet_jacobi_sum_gauss(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)

.. function:: void acb_dirichlet_jacobi_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1,  const dirichlet_char_t chi2, slong prec)

.. function:: void acb_dirichlet_jacobi_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, ulong b, slong prec)

   Computes the Jacobi sum

   .. math::

      J_q(a,b) = \sum_{x \bmod q} \chi_q(a, x)\chi_q(b, 1-x)

   - the *naive* version computes the sum as defined.

   - the *factor* version writes it as a product of local Jacobi sums

   - the *gauss* version assumes `ab` is primitive and uses the formula
     `J_q(a,b)G_q(ab) = G_q(a)G_q(b)`

   - the default version automatically combines the above methods.

   - the *ui* version only takes the Conrey numbers *a* and *b* as parameters.

Theta sums
-------------------------------------------------------------------------------

We call *theta series* of a Dirichlet character the quadratic series

.. math::

   \Theta_q(a) = \sum_{n\geq 0} \chi_q(a, n) n^p x^{n^2}

where `p` is the parity of the character `\chi_q(a,\cdot)`.

For `\Re(t)>0` we write `x(t)=\exp(-\frac{\pi}{N}t^2)` and define

.. math::

   \Theta_q(a,t) = \sum_{n\geq 0} \chi_q(a, n) x(t)^{n^2}.

.. function:: void acb_dirichlet_chi_theta_arb(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t t, slong prec)

.. function:: void acb_dirichlet_ui_theta_arb(acb_t res, const dirichlet_group_t G, ulong a, const arb_t t, slong prec)

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


.. function:: void acb_dirichlet_dft_conrey(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

   Compute the DFT of *v* using Conrey indices.
   This function assumes *v* and *w* are vectors
   of size *G->phi_q*, whose values correspond to a lexicographic ordering
   of Conrey logs (as obtained using :func:`dirichlet_conrey_next` or
   by :func:`dirichlet_index_conrey`).

   For example, if `q=15`, the Conrey elements are stored in following
   order

   =======  =============  =====================
    index    log = [e,f]     number = 7^e 11^f
   =======  =============  =====================
      0       [0, 0]        1
      1       [0, 1]        7
      2       [0, 2]        4
      3       [0, 3]        13
      4       [0, 4]        1
      5       [1, 0]        11
      6       [1, 1]        2
      7       [1, 2]        14
      8       [1, 3]        8
      9       [1, 4]        11
   =======  =============  =====================

.. function:: void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

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

.. function:: void acb_dirichlet_root_number_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

.. function:: void acb_dirichlet_root_number(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

   Sets *res* to the root number `\epsilon(\chi)` for a primitive character *chi*,
   which appears in the functional equation (where `p` is the parity of `\chi`):

   .. math::

      \left(\frac{q}{\pi}\right)^{\frac{s+p}2}\Gamma\left(\frac{s+p}2\right) L(s, \chi) = \epsilon(\chi) \left(\frac{q}{\pi}\right)^{\frac{1-s+p}2}\Gamma\left(\frac{1-s+p}2\right) L(1 - s, \overline\chi)

   - The *theta* variant uses the evaluation at `t=1` of the Theta series.

   - The default version computes it via the gauss sum.

.. function:: void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    Computes `L(s,\chi)` using decomposition in terms of the Hurwitz zeta function

    .. math::

        L(s,\chi) = q^{-s}\sum_{k=1}^{q-1} \chi(k) \,\zeta\!\left(s,\frac kq\right).

    If `s = 1` and `\chi` is non-principal, the deflated Hurwitz zeta function
    is used to avoid poles.

    This formula is slow for large *q*.

.. function:: void acb_dirichlet_l_euler_product(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    Computes `L(s,\chi)` directly using the Euler product. This is
    efficient if *s* has large positive real part. As implemented, this
    function only gives a finite result if `\operatorname{re}(s) \ge 2`.

.. function:: void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    Computes `L(s,\chi)` using a default choice of algorithm.

.. function:: void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const dirichlet_group_t G, slong prec)

    Compute all values `L(s,\chi)` for `\chi` mod `q`, by Hurwitz formula and
    discrete Fourier transform.
    *res* is assumed to have length *G->phi_q* and values are stored by lexicographically ordered
    Conrey logs. See :func:`acb_dirichlet_dft_conrey`.
