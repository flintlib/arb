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

