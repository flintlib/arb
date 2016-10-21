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

Truncated L-series and power sums
-------------------------------------------------------------------------------

.. function:: void acb_dirichlet_powsum_term(acb_ptr res, arb_t log_prev, ulong * prev, const acb_t s, ulong k, int integer, int critical_line, slong len, slong prec)

    Sets *res* to `k^{-(s+x)}` as a power series in *x* truncated to length *len*.
    The flags *integer* and *critical_line* respectively specify optimizing
    for *s* being an integer or having real part 1/2.

    On input *log_prev* should contain the natural logarithm of the integer
    at *prev*. If *prev* is close to *k*, this can be used to speed up
    computations. If `\log(k)` is computed internally by this function, then
    *log_prev* is overwritten by this value, and the integer at *prev* is
    overwritten by *k*, allowing *log_prev* to be recycled for the next
    term when evaluating a power sum.

.. function:: void acb_dirichlet_powsum_sieved(acb_ptr res, const acb_t s, ulong n, slong len, slong prec)

    Sets *res* to `\sum_{k=1}^n k^{-(s+x)}`
    as a power series in *x* truncated to length *len*.
    This function stores a table of powers that have already been calculated,
    computing `(ij)^r` as `i^r j^r` whenever `k = ij` is
    composite. As a further optimization, it groups all even `k` and
    evaluates the sum as a polynomial in `2^{-(s+x)}`.
    This scheme requires about `n / \log n` powers, `n / 2` multiplications,
    and temporary storage of `n / 6` power series. Due to the extra
    power series multiplications, it is only faster than the naive
    algorithm when *len* is small.

.. function:: void acb_dirichlet_powsum_smooth(acb_ptr res, const acb_t s, ulong n, slong len, slong prec)

    Sets *res* to `\sum_{k=1}^n k^{-(s+x)}`
    as a power series in *x* truncated to length *len*.
    This function performs partial sieving by adding multiples of 5-smooth *k*
    into separate buckets. Asymptotically, this requires computing 4/15
    of the powers, which is slower than *sieved*, but only requires
    logarithmic extra space. It is also faster for large *len*, since most
    power series multiplications are traded for additions.
    A slightly bigger gain for larger *n* could be achieved by using more
    small prime factors, at the expense of space.

Hurwitz zeta function
-------------------------------------------------------------------------------

.. type:: acb_dirichlet_hurwitz_precomp_struct

.. type:: acb_dirichlet_hurwitz_precomp_t

.. function:: void acb_dirichlet_hurwitz_precomp_init(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, ulong A, ulong K, ulong N, slong prec)

    Precomputes a grid of Taylor polynomials for fast evaluation of
    `\zeta(s,a)` on `a \in (0,1]` with fixed *s*.
    *A* is the initial shift to apply to *a*, *K* is the number of Taylor terms,
    *N* is the number of grid points.  The precomputation requires *NK*
    evaluations of the Hurwitz zeta function, and each subsequent evaluation
    requires *2K* simple arithmetic operations (polynomial evaluation) plus
    *A* powers. As *K* grows, the error is at most `O(1/(2AN)^K)`.

    We require that *A*, *K* and *N* are all positive. Moreover, for a finite
    error bound, we require `K+\operatorname{re}(s) > 1`.
    To avoid an initial "bump" that steals precision
    and slows convergence, *AN* should be at least roughly as large as `|s|`,
    e.g. it is a good idea to have at least `AN > 0.5 |s|`.

.. function:: void acb_dirichlet_hurwitz_precomp_clear(acb_dirichlet_hurwitz_precomp_t pre)

    Clears the precomputed data.

.. function:: void acb_dirichlet_hurwitz_precomp_bound(mag_t res, const acb_t s, ulong A, ulong K, ulong N)

    Computes an upper bound for the truncation error (not accounting for
    roundoff error) when evaluating `\zeta(s,a)` with precomputation parameters
    *A*, *K*, *N*, assuming that `0 < a \le 1`.
    For details, see :ref:`algorithms_hurwitz`.

.. function:: void acb_dirichlet_hurwitz_precomp_eval(acb_t res, const acb_dirichlet_hurwitz_precomp_t pre, ulong p, ulong q, slong prec)

    Evaluates `\zeta(s,p/q)` using precomputed data, assuming that `0 < p/q \le 1`.

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

        L(s,\chi) = q^{-s}\sum_{k=1}^q \chi(k) \,\zeta\!\left(s,\frac kq\right).

    If `s = 1` and `\chi` is non-principal, the deflated Hurwitz zeta function
    is used to avoid poles.

    This formula is slow for large *q*.

.. function:: void acb_dirichlet_l_euler_product(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

.. function:: void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s, const signed char * chi, int mod, int reciprocal, slong prec)

    Computes `L(s,\chi)` directly using the Euler product. This is
    efficient if *s* has large positive real part. As implemented, this
    function only gives a finite result if `\operatorname{re}(s) \ge 2`.

    An error bound is computed via :func:`mag_hurwitz_zeta_uiui`.
    If *s* is complex, replace it with its real part. Since

    .. math ::

        \frac{1}{L(s,\chi)} = \prod_{p} \left(1 - \frac{\chi(p)}{p^s}\right)
                = \sum_{k=1}^{\infty} \frac{\mu(k)\chi(k)}{k^s}

    and the truncated product gives all smooth-index terms in the series, we have

    .. math ::

        \left|\prod_{p < N} \left(1 - \frac{\chi(p)}{p^s}\right) - \frac{1}{L(s,\chi)}\right|
        \le \sum_{k=N}^{\infty} \frac{1}{k^s} = \zeta(s,N).

    The underscore version specialized for integer *s* assumes that `\chi` is
    a real Dirichlet character given by the explicit list *chi* of character
    values at 0, 1, ..., *mod* - 1. If *reciprocal* is set, it computes
    `1 / L(s,\chi)` (this is faster if the reciprocal can be used directly).

.. function:: void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    Computes `L(s,\chi)` using a default choice of algorithm.

.. function:: void acb_dirichlet_l_jet(acb_ptr res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)

    Computes the Taylor expansion of `L(s,\chi)` to length *len*,
    i.e. `L(s), L'(s), \ldots, L^{(len-1)}(s) / (len-1)!`.
    If *deflate* is set, computes the expansion of

    .. math ::

        L(s,\chi) - \frac{\sum_{k=1}^q \chi(k)}{(s-1)q}

    instead. If *chi* is a principal character, then this has the effect of
    subtracting the pole with residue `\sum_{k=1}^q \chi(k) = \phi(q) / q`
    that is located at `s = 1`. In particular, when evaluated at `s = 1`, this
    gives the regular part of the Laurent expansion.
    When *chi* is non-principal, *deflate* has no effect.

