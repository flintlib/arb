.. _acb-hypgeom:

**acb_hypgeom.h** -- hypergeometric functions of complex variables
==================================================================================

The generalized hypergeometric function is formally defined by

.. math ::

    {}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z) =
    \sum_{k=0}^\infty \frac{(a_1)_k\dots(a_p)_k}{(b_1)_k\dots(b_q)_k} \frac {z^k} {k!}.

It can be interpreted using analytic continuation or regularization
when the sum does not converge.
In a looser sense, we understand "hypergeometric functions" to be
linear combinations of generalized hypergeometric functions
with prefactors that are products of exponentials, powers, and gamma functions.

Convergent series
-------------------------------------------------------------------------------

In this section, we define

.. math ::

    T(k) = \frac{\prod_{i=0}^{p-1} (a_i)_k}{\prod_{i=0}^{q-1} (b_i)_k} z^k

and

.. math ::

    {}_pf_{q}(a_0,\ldots,a_{p-1}; b_0 \ldots b_{q-1}; z) = {}_{p+1}F_{q}(a_0,\ldots,a_{p-1},1; b_0 \ldots b_{q-1}; z) = \sum_{k=0}^{\infty} T(k)

For the conventional generalized hypergeometric function
`{}_pF_{q}`, compute  `{}_pf_{q+1}` with the explicit parameter `b_q = 1`,
or remove a 1 from the `a_i` parameters if there is one.

.. function:: void acb_hypgeom_pfq_bound_factor(mag_t C, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, ulong n)

    Computes a factor *C* such that
    `\left|\sum_{k=n}^{\infty} T(k)\right| \le C |T(n)|`.
    See :ref:`algorithms_hypergeometric_convergent`.
    As currently implemented, the bound becomes infinite when `n` is
    too small, even if the series converges.

.. function:: slong acb_hypgeom_pfq_choose_n(acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec)

    Heuristically attempts to choose a number of terms *n* to
    sum of a hypergeometric series at a working precision of *prec* bits.

    Uses double precision arithmetic internally. As currently implemented,
    it can fail to produce a good result if the parameters are extremely
    large or extremely close to nonpositive integers.

    Numerical cancellation is assumed to be significant, so truncation
    is done when the current term is *prec* bits
    smaller than the largest encountered term.

    This function will also attempt to pick a reasonable
    truncation point for divergent series.

.. function:: void acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

.. function:: void acb_hypgeom_pfq_sum_rs(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

.. function:: void acb_hypgeom_pfq_sum_bs(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

.. function:: void acb_hypgeom_pfq_sum_fme(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

.. function:: void acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

    Computes `s = \sum_{k=0}^{n-1} T(k)` and `t = T(n)`.
    Does not allow aliasing between input and output variables.
    We require `n \ge 0`.

    The *forward* version computes the sum using forward
    recurrence.

    The *bs* version computes the sum using binary splitting.

    The *rs* version computes the sum in reverse order
    using rectangular splitting. It only computes a
    magnitude bound for the value of *t*.

    The *fme* version uses fast multipoint evaluation.

    The default version automatically chooses an algorithm
    depending on the inputs.

.. function:: void acb_hypgeom_pfq_sum_bs_invz(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t w, slong n, slong prec)

.. function:: void acb_hypgeom_pfq_sum_invz(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, const acb_t w, slong n, slong prec)

    Like :func:`acb_hypgeom_pfq_sum`, but taking advantage of
    `w = 1/z` possibly having few bits.

.. function:: void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

    Computes

    .. math ::

        {}_pf_{q}(z)
            = \sum_{k=0}^{\infty} T(k)
            = \sum_{k=0}^{n-1} T(k) + \varepsilon

    directly from the defining series, including a rigorous bound for
    the truncation error `\varepsilon` in the output.

    If  `n < 0`, this function chooses a number of terms automatically
    using :func:`acb_hypgeom_pfq_choose_n`.

.. function:: void acb_hypgeom_pfq_series_direct(acb_poly_t res, const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, int regularized, slong n, slong len, slong prec)

    Computes `{}_pf_{q}(z)` directly using the defining series, given
    parameters and argument that are power series.
    The result is a power series of length *len*.

    An error bound is computed automatically as a function of the number
    of terms *n*. If `n < 0`, the number of terms is chosen
    automatically.

    If *regularized* is set, the regularized hypergeometric function
    is computed instead.

Asymptotic series
-------------------------------------------------------------------------------

`U(a,b,z)` is the confluent hypergeometric function of the second
kind with the principal branch cut, and `U^{*} = z^a U(a,b,z)`.
For details about how error bounds are computed,
see :ref:`algorithms_hypergeometric_asymptotic_confluent`.

.. function:: void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong n, slong prec)

    Sets *res* to `U^{*}(a,b,z)` computed using *n* terms of the asymptotic series,
    with a rigorous bound for the error included in the output.
    We require `n \ge 0`.

.. function:: int acb_hypgeom_u_use_asymp(const acb_t z, slong prec)

    Heuristically determines whether the asymptotic series can be used
    to evaluate `U(a,b,z)` to *prec* accurate bits (assuming that *a* and *b*
    are small).

Generalized hypergeometric function
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_pfq(acb_poly_t res, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, int regularized, slong prec)

    Computes the generalized hypergeometric function `{}_pF_{q}(z)`,
    or the regularized version if *regularized* is set.

    This function automatically delegates to a specialized implementation
    when the order (*p*, *q*) is one of (0,0), (1,0), (0,1), (1,1), (2,1).
    Otherwise, it falls back to direct summation.

    While this is a top-level function meant to take care of special cases
    automatically, it does not generally perform the optimization
    of deleting parameters that appear in both *a* and *b*. This can be
    done ahead of time by the user in applications where duplicate
    parameters are likely to occur.

Confluent hypergeometric functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_u_1f1_series(acb_poly_t res, const acb_poly_t a, const acb_poly_t b, const acb_poly_t z, slong len, slong prec)

    Computes `U(a,b,z)` as a power series truncated to length *len*,
    given `a, b, z \in \mathbb{C}[[x]]`.
    If `b[0] \in \mathbb{Z}`, it computes one extra derivative and removes
    the singularity (it is then assumed that `b[1] \ne 0`).
    As currently implemented, the output is indeterminate if `b` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)

    Computes `U(a,b,z)` as a sum of two convergent hypergeometric series.
    If `b \in \mathbb{Z}`, it computes
    the limit value via :func:`acb_hypgeom_u_1f1_series`.
    As currently implemented, the output is indeterminate if `b` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)

    Computes `U(a,b,z)` using an automatic algorithm choice. The
    function :func:`acb_hypgeom_u_asymp` is used
    if `a` or `a-b+1` is a nonpositive integer (in which
    case the asymptotic series terminates), or if *z* is sufficiently large.
    Otherwise :func:`acb_hypgeom_u_1f1` is used.

.. function:: void acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)

    Computes the confluent hypergeometric function
    `M(a,b,z) = {}_1F_1(a,b,z)`, or
    `\mathbf{M}(a,b,z) = \frac{1}{\Gamma(b)} {}_1F_1(a,b,z)` if *regularized*
    is set.

.. function:: void acb_hypgeom_0f1_asymp(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_0f1_direct(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec)

    Computes the confluent hypergeometric function
    `{}_0F_1(a,z)`, or `\frac{1}{\Gamma(a)} {}_0F_1(a,z)` if *regularized*
    is set, using asymptotic expansions, direct summation,
    or an automatic algorithm choice.
    The *asymp* version uses the asymptotic expansions of Bessel
    functions, together with the connection formulas

    .. math ::

        \frac{{}_0F_1(a,z)}{\Gamma(a)} = (-z)^{(1-a)/2} J_{a-1}(2 \sqrt{-z}) =
                                         z^{(1-a)/2} I_{a-1}(2 \sqrt{z}).

    The Bessel-*J* function is used in the left half-plane and the
    Bessel-*I* function is used in the right half-plane, to avoid loss
    of accuracy due to evaluating the square root on the branch cut.

Error functions and Fresnel integrals
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_erf_propagated_error(mag_t re, mag_t im, const acb_t z)

    Sets *re* and *im* to upper bounds for the error in the real and imaginary
    part resulting from approximating the error function of *z* by
    the error function evaluated at the midpoint of *z*. Uses
    the first derivative.

.. function:: void acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_erf_asymp(acb_t res, const acb_t z, int complementary, slong prec, slong prec2)

    Computes the error function respectively using

    .. math ::

        \operatorname{erf}(z) &= \frac{2z}{\sqrt{\pi}}
            {}_1F_1(\tfrac{1}{2}, \tfrac{3}{2}, -z^2)

        \operatorname{erf}(z) &= \frac{2z e^{-z^2}}{\sqrt{\pi}}
            {}_1F_1(1, \tfrac{3}{2}, z^2)

        \operatorname{erf}(z) &= \frac{z}{\sqrt{z^2}}
            \left(1 - \frac{e^{-z^2}}{\sqrt{\pi}}
            U(\tfrac{1}{2}, \tfrac{1}{2}, z^2)\right) =
            \frac{z}{\sqrt{z^2}} - \frac{e^{-z^2}}{z \sqrt{\pi}}
            U^{*}(\tfrac{1}{2}, \tfrac{1}{2}, z^2).

    The *asymp* version takes a second precision to use for the *U* term.
    It also takes an extra flag *complementary*, computing the complementary
    error function if set.

.. function:: void acb_hypgeom_erf(acb_t res, const acb_t z, slong prec)

    Computes the error function using an automatic algorithm choice.
    If *z* is too small to use the asymptotic expansion, a working precision
    sufficient to circumvent cancellation in the hypergeometric series is
    determined automatically, and a bound for the propagated error is
    computed with :func:`acb_hypgeom_erf_propagated_error`.

.. function:: void _acb_hypgeom_erf_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_erf_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the error function of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_erfc(acb_t res, const acb_t z, slong prec)

    Computes the complementary error function
    `\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)`.
    This function avoids catastrophic cancellation for large positive *z*.

.. function:: void _acb_hypgeom_erfc_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_erfc_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the complementary error function of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_erfi(acb_t res, const acb_t z, slong prec)

    Computes the imaginary error function
    `\operatorname{erfi}(z) = -i\operatorname{erf}(iz)`. This is a trivial wrapper
    of :func:`acb_hypgeom_erf`.

.. function:: void _acb_hypgeom_erfi_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_erfi_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the imaginary error function of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_fresnel(acb_t res1, acb_t res2, const acb_t z, int normalized, slong prec)

    Sets *res1* to the Fresnel sine integral `S(z)` and *res2* to
    the Fresnel cosine integral `C(z)`. Optionally, just a single function
    can be computed by passing *NULL* as the other output variable.
    The definition `S(z) = \int_0^z \sin(t^2) dt` is used if *normalized* is 0,
    and `S(z) = \int_0^z \sin(\tfrac{1}{2} \pi t^2) dt` is used if
    *normalized* is 1 (the latter is the Abramowitz & Stegun convention).
    `C(z)` is defined analogously.

.. function:: void _acb_hypgeom_fresnel_series(acb_ptr res1, acb_ptr res2, acb_srcptr z, slong zlen, int normalized, slong len, slong prec)

.. function:: void acb_hypgeom_fresnel_series(acb_poly_t res1, acb_poly_t res2, const acb_poly_t z, int normalized, slong len, slong prec)

    Sets *res1* to the Fresnel sine integral and *res2* to the Fresnel
    cosine integral of the power series *z*, truncated to length *len*.
    Optionally, just a single function can be computed by passing *NULL*
    as the other output variable.

Bessel functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the Bessel function of the first kind
    via :func:`acb_hypgeom_u_asymp`.
    For all complex `\nu, z`, we have

    .. math ::

        J_{\nu}(z) = \frac{z^{\nu}}{2^{\nu} e^{iz} \Gamma(\nu+1)}
            {}_1F_1(\nu+\tfrac{1}{2}, 2\nu+1, 2iz) = A_{+} B_{+} + A_{-} B_{-}

    where

    .. math ::

        A_{\pm} = z^{\nu} (z^2)^{-\tfrac{1}{2}-\nu} (\mp i z)^{\tfrac{1}{2}+\nu} (2 \pi)^{-1/2} = (\pm iz)^{-1/2-\nu} z^{\nu} (2 \pi)^{-1/2}

    .. math ::

        B_{\pm} = e^{\pm i z} U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, \mp 2iz).

    Nicer representations of the factors `A_{\pm}` can be given depending conditionally
    on the parameters. If `\nu + \tfrac{1}{2} = n \in \mathbb{Z}`, we have
    `A_{\pm} = (\pm i)^{n} (2 \pi z)^{-1/2}`.
    And if `\operatorname{Re}(z) > 0`, we have `A_{\pm} = \exp(\mp i [(2\nu+1)/4] \pi) (2 \pi z)^{-1/2}`.

.. function:: void acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the Bessel function of the first kind from

    .. math ::

        J_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
                     {}_0F_1\left(\nu+1, -\frac{z^2}{4}\right).

.. function:: void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the Bessel function of the first kind `J_{\nu}(z)` using
    an automatic algorithm choice.

.. function:: void acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the Bessel function of the second kind `Y_{\nu}(z)` from the
    formula

    .. math ::

        Y_{\nu}(z) = \frac{\cos(\nu \pi) J_{\nu}(z) - J_{-\nu}(z)}{\sin(\nu \pi)}

    unless `\nu = n` is an integer in which case the limit value

    .. math ::

        Y_n(z) = -\frac{2}{\pi} \left( i^n K_n(iz) +
            \left[\log(iz)-\log(z)\right] J_n(z) \right)

    is computed.
    As currently implemented, the output is indeterminate if `\nu` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_bessel_jy(acb_t res1, acb_t res2, const acb_t nu, const acb_t z, slong prec)

    Sets *res1* to `J_{\nu}(z)` and *res2* to `Y_{\nu}(z)`, computed
    simultaneously. From these values, the user can easily
    construct the Bessel functions of the third kind (Hankel functions)
    `H_{\nu}^{(1)}(z), H_{\nu}^{(2)}(z) = J_{\nu}(z) \pm i Y_{\nu}(z)`.

.. function:: void acb_hypgeom_bessel_i_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)

.. function:: void acb_hypgeom_bessel_i_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec)

.. function:: void acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the modified Bessel function of the first kind
    `I_{\nu}(z) = z^{\nu} (iz)^{-\nu} J_{\nu}(iz)` respectively using
    asymptotic series (see :func:`acb_hypgeom_bessel_j_asymp`),
    the convergent series

    .. math ::

        I_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
                     {}_0F_1\left(\nu+1, \frac{z^2}{4}\right),

    or an automatic algorithm choice.

.. function:: void acb_hypgeom_bessel_k_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the modified Bessel function of the second kind via
    via :func:`acb_hypgeom_u_asymp`. For all `\nu` and all `z \ne 0`, we have

    .. math ::

        K_{\nu}(z) = \left(\frac{2z}{\pi}\right)^{-1/2} e^{-z}
            U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, 2z).

.. function:: void acb_hypgeom_bessel_k_0f1_series(acb_poly_t res, const acb_poly_t nu, const acb_poly_t z, slong len, slong prec)

    Computes the modified Bessel function of the second kind `K_{\nu}(z)`
    as a power series truncated to length *len*,
    given `\nu, z \in \mathbb{C}[[x]]`. Uses the formula

    .. math ::

        K_{\nu}(z) = \frac{1}{2} \frac{\pi}{\sin(\pi \nu)} \left[
                    \left(\frac{z}{2}\right)^{-\nu}
                        {}_0{\widetilde F}_1\left(1-\nu, \frac{z^2}{4}\right)
                     -
                     \left(\frac{z}{2}\right)^{\nu}
                         {}_0{\widetilde F}_1\left(1+\nu, \frac{z^2}{4}\right)
                    \right].

    If `\nu[0] \in \mathbb{Z}`, it computes one extra derivative and removes
    the singularity (it is then assumed that `\nu[1] \ne 0`).
    As currently implemented, the output is indeterminate if `\nu[0]` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_bessel_k_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the modified Bessel function of the second kind from

    .. math ::

        K_{\nu}(z) = \frac{1}{2} \left[
                    \left(\frac{z}{2}\right)^{-\nu}
                        \Gamma(\nu)
                        {}_0F_1\left(1-\nu, \frac{z^2}{4}\right)
                     -
                     \left(\frac{z}{2}\right)^{\nu}
                         \frac{\pi}{\nu \sin(\pi \nu) \Gamma(\nu)}
                         {}_0F_1\left(\nu+1, \frac{z^2}{4}\right)
                    \right]

    if `\nu \notin \mathbb{Z}`. If `\nu \in \mathbb{Z}`, it computes
    the limit value via :func:`acb_hypgeom_bessel_k_0f1_series`.
    As currently implemented, the output is indeterminate if `\nu` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, slong prec)

    Computes the modified Bessel function of the second kind `K_{\nu}(z)` using
    an automatic algorithm choice.

Airy functions
-------------------------------------------------------------------------------

The Airy functions are linearly independent solutions of the
differential equation `y'' - zy = 0`. All solutions are entire functions.
The standard solutions are denoted `\operatorname{Ai}(z), \operatorname{Bi}(z)`.
For negative *z*, both functions are oscillatory. For positive *z*,
the first function decreases exponentially while the second increases
exponentially.

The Airy functions can be expressed in terms of Bessel functions of fractional
order, but this is inconvenient since such formulas
only hold piecewise (due to the Stokes phenomenon). Computation of the
Airy functions can also be optimized more than Bessel functions in general.
We therefore provide a dedicated interface for evaluating Airy functions.

The following methods optionally compute
`(\operatorname{Ai}(z), \operatorname{Ai}'(z), \operatorname{Bi}(z), \operatorname{Bi}'(z))`
simultaneously. Any of the four function values can be omitted by passing
*NULL* for the unwanted output variables, speeding up the evaluation.

.. function:: void acb_hypgeom_airy_direct(acb_t ai, acb_t ai_prime, acb_t bi, acb_t bi_prime, const acb_t z, slong n, slong prec)

    Computes the Airy functions using direct series expansions truncated at *n* terms.
    Error bounds are included in the output.

.. function:: void acb_hypgeom_airy_asymp(acb_t ai, acb_t ai_prime, acb_t bi, acb_t bi_prime, const acb_t z, slong n, slong prec)

    Computes the Airy functions using asymptotic expansions truncated at *n* terms.
    Error bounds are included in the output.
    For details about how the error bounds are computed, see
    :ref:`algorithms_hypergeometric_asymptotic_airy`.

.. function:: void acb_hypgeom_airy_bound(mag_t ai, mag_t ai_prime, mag_t bi, mag_t bi_prime, const acb_t z)

    Computes bounds for the Airy functions using first-order asymptotic
    expansions together with error bounds. This function uses some
    shortcuts to make it slightly faster than calling
    :func:`acb_hypgeom_airy_asymp` with `n = 1`.

.. function:: void acb_hypgeom_airy(acb_t ai, acb_t ai_prime, acb_t bi, acb_t bi_prime, const acb_t z, slong prec)

    Computes Airy functions using an automatic algorithm choice.

    We use :func:`acb_hypgeom_airy_asymp` whenever this gives full accuracy
    and :func:`acb_hypgeom_airy_direct` otherwise.
    In the latter case, we first use hardware double precision arithmetic to
    determine an accurate estimate of the working precision needed
    to compute the Airy functions accurately for given *z*. This estimate is
    obtained by comparing the leading-order asymptotic estimate of the Airy
    functions with the magnitude of the largest term in the power series.
    The estimate is generic in the sense that it does not take into account
    vanishing near the roots of the functions.
    We subsequently evaluate the power series at the midpoint of *z* and
    bound the propagated error using derivatives. Derivatives are
    bounded using :func:`acb_hypgeom_airy_bound`.

.. function:: void acb_hypgeom_airy_jet(acb_ptr ai, acb_ptr bi, const acb_t z, slong len, slong prec)

    Writes to *ai* and *bi* the respective Taylor expansions of the Airy functions
    at the point *z*, truncated to length *len*.
    Either of the outputs can be *NULL* to avoid computing that function.
    The variable *z* is not allowed to be aliased with the outputs.
    To simplify the implementation, this method does not compute the
    series expansions of the primed versions directly; these are
    easily obtained by computing one extra coefficient and differentiating
    the output with :func:`_acb_poly_derivative`.

.. function:: void _acb_hypgeom_airy_series(acb_ptr ai, acb_ptr ai_prime, acb_ptr bi, acb_ptr bi_prime, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_airy_series(acb_poly_t ai, acb_poly_t ai_prime, acb_poly_t bi, acb_poly_t bi_prime, const acb_poly_t z, slong len, slong prec)

    Computes the Airy functions evaluated at the power series *z*,
    truncated to length *len*. As with the other Airy methods, any of the
    outputs can be *NULL*.

Incomplete gamma functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_gamma_upper_singular(acb_t res, slong s, const acb_t z, int regularized, slong prec)

.. function:: void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

    If *regularized* is 0, computes the upper incomplete gamma function
    `\Gamma(s,z)`.

    If *regularized* is 1, computes the regularized upper incomplete
    gamma function `Q(s,z) = \Gamma(s,z) / \Gamma(s)`.

    If *regularized* is 2, computes the generalized exponential integral
    `z^{-s} \Gamma(s,z) = E_{1-s}(z)` instead (this option is mainly
    intended for internal use; :func:`acb_hypgeom_expint` is the intended
    interface for computing the exponential integral).

    The different methods respectively implement the formulas

    .. math ::

        \Gamma(s,z) = e^{-z} U(1-s,1-s,z)

    .. math ::

        \Gamma(s,z) = \Gamma(s) - \frac{z^s}{s} {}_1F_1(s, s+1, -z)

    .. math ::

        \Gamma(s,z) = \Gamma(s) - \frac{z^s e^{-z}}{s} {}_1F_1(1, s+1, z)

    .. math ::

        \Gamma(s,z) = \frac{(-1)^n}{n!} (\psi(n+1) - \log(z))
                    + \frac{(-1)^n}{(n+1)!} z \, {}_2F_2(1,1,2,2+n,-z)
                    - z^{-n} \sum_{k=0}^{n-1} \frac{(-z)^k}{(k-n) k!},
                    \quad n = -s \in \mathbb{Z}_{\ge 0}

    and an automatic algorithm choice. The automatic version also handles
    other special input such as `z = 0` and `s = 1, 2, 3`.
    The *singular* version evaluates the finite sum directly and therefore
    assumes that *s* is not too large.

.. function:: void _acb_hypgeom_gamma_upper_series(acb_ptr res, acb_t s, acb_srcptr z, slong zlen, int regularized, slong n, slong prec)

.. function:: void acb_hypgeom_gamma_upper_series(acb_poly_t res, const acb_t s, const acb_poly_t z, int regularized, slong n, slong prec)

    Sets *res* to an upper incomplete gamma function where *s* is
    a constant and *z* is a power series, truncated to length *n*.
    The *regularized* argument has the same interpretation as in
    :func:`acb_hypgeom_gamma_upper`.


.. function:: void acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

    If *regularized* is 0, computes the lower incomplete gamma function
    `\gamma(s,z) = \frac{z^s}{s} {}_1F_1(s, s+1, -z)`.

    If *regularized* is 1, computes the regularized lower incomplete
    gamma function `P(s,z) = \gamma(s,z) / \Gamma(s)`.

    If *regularized* is 2, computes a further regularized lower incomplete
    gamma function `\gamma^{*}(s,z) = z^{-s} P(s,z)`.

.. function:: void _acb_hypgeom_gamma_lower_series(acb_ptr res, acb_t s, acb_srcptr z, slong zlen, int regularized, slong n, slong prec)

.. function:: void acb_hypgeom_gamma_lower_series(acb_poly_t res, const acb_t s, const acb_poly_t z, int regularized, slong n, slong prec)

    Sets *res* to an lower incomplete gamma function where *s* is
    a constant and *z* is a power series, truncated to length *n*.
    The *regularized* argument has the same interpretation as in
    :func:`acb_hypgeom_gamma_lower`.

Exponential and trigonometric integrals
-------------------------------------------------------------------------------

The branch cut conventions of the following functions match Mathematica.

.. function:: void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, slong prec)

    Computes the generalized exponential integral `E_s(z)`. This is a
    trivial wrapper of :func:`acb_hypgeom_gamma_upper`.

.. function:: void acb_hypgeom_ei_asymp(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_ei_2f2(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_ei(acb_t res, const acb_t z, slong prec)

    Computes the exponential integral `\operatorname{Ei}(z)`, respectively
    using

    .. math ::

        \operatorname{Ei}(z) = -e^z U(1,1,-z) - \log(-z)
            + \frac{1}{2} \left(\log(z) - \log\left(\frac{1}{z}\right) \right)

    .. math ::

        \operatorname{Ei}(z) = z {}_2F_2(1, 1; 2, 2; z) + \gamma
            + \frac{1}{2} \left(\log(z) - \log\left(\frac{1}{z}\right) \right)

    and an automatic algorithm choice.

.. function:: void _acb_hypgeom_ei_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_ei_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the exponential integral of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_si_asymp(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_si_1f2(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_si(acb_t res, const acb_t z, slong prec)

    Computes the sine integral `\operatorname{Si}(z)`, respectively
    using

    .. math ::

        \operatorname{Si}(z) = \frac{i}{2} \left[
            e^{iz} U(1,1,-iz) - e^{-iz} U(1,1,iz) + 
            \log(-iz) - \log(iz) \right]

    .. math ::

        \operatorname{Si}(z) = z {}_1F_2(\tfrac{1}{2}; \tfrac{3}{2}, \tfrac{3}{2}; -\tfrac{z^2}{4})

    and an automatic algorithm choice.

.. function:: void _acb_hypgeom_si_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_si_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the sine integral of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_ci_asymp(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_ci_2f3(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_ci(acb_t res, const acb_t z, slong prec)

    Computes the cosine integral `\operatorname{Ci}(z)`, respectively
    using

    .. math ::

        \operatorname{Ci}(z) = \log(z) - \frac{1}{2} \left[
            e^{iz} U(1,1,-iz) + e^{-iz} U(1,1,iz) + 
            \log(-iz) + \log(iz) \right]

    .. math ::

        \operatorname{Ci}(z) = -\tfrac{z^2}{4}
            {}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; -\tfrac{z^2}{4})
            + \log(z) + \gamma

    and an automatic algorithm choice.

.. function:: void _acb_hypgeom_ci_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_ci_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the cosine integral of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_shi(acb_t res, const acb_t z, slong prec)

    Computes the hyperbolic sine integral
    `\operatorname{Shi}(z) = -i \operatorname{Si}(iz)`.
    This is a trivial wrapper of :func:`acb_hypgeom_si`.

.. function:: void _acb_hypgeom_shi_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_shi_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the hyperbolic sine integral of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_chi_asymp(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_chi_2f3(acb_t res, const acb_t z, slong prec)

.. function:: void acb_hypgeom_chi(acb_t res, const acb_t z, slong prec)

    Computes the hyperbolic cosine integral `\operatorname{Chi}(z)`, respectively
    using

    .. math ::

        \operatorname{Chi}(z) = -\frac{1}{2} \left[
            e^{z} U(1,1,-z) + e^{-z} U(1,1,z) + 
            \log(-z) - \log(z) \right]

    .. math ::

        \operatorname{Chi}(z) = \tfrac{z^2}{4}
            {}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; \tfrac{z^2}{4})
            + \log(z) + \gamma

    and an automatic algorithm choice.

.. function:: void _acb_hypgeom_chi_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void acb_hypgeom_chi_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

    Computes the hyperbolic cosine integral of the power series *z*,
    truncated to length *len*.

.. function:: void acb_hypgeom_li(acb_t res, const acb_t z, int offset, slong prec)

    If *offset* is zero, computes the logarithmic integral
    `\operatorname{li}(z) = \operatorname{Ei}(\log(z))`.

    If *offset* is nonzero, computes the offset logarithmic integral
    `\operatorname{Li}(z) = \operatorname{li}(z) - \operatorname{li}(2)`.

.. function:: void _acb_hypgeom_li_series(acb_ptr res, acb_srcptr z, slong zlen, int offset, slong len, slong prec)

.. function:: void acb_hypgeom_li_series(acb_poly_t res, const acb_poly_t z, int offset, slong len, slong prec)

    Computes the logarithmic integral (optionally the offset version)
    of the power series *z*, truncated to length *len*.

Gauss hypergeometric function
-------------------------------------------------------------------------------

The following methods compute the Gauss hypergeometric function

.. math ::

    F(z) = {}_2F_1(a,b,c,z) = \sum_{k=0}^{\infty} \frac{(a)_k (b)_k}{(c)_k} \frac{z^k}{k!}

or the regularized version
`\operatorname{\mathbf{F}}(z) = \operatorname{\mathbf{F}}(a,b,c,z) = {}_2F_1(a,b,c,z) / \Gamma(c)`
if the flag *regularized* is set.

.. function:: void acb_hypgeom_2f1_continuation(acb_t res0, acb_t res1, const acb_t a, const acb_t b, const acb_t c, const acb_t z0, const acb_t z1, const acb_t f0, const acb_t f1, slong prec)

    Given `F(z_0), F'(z_0)` in *f0*, *f1*, sets *res0* and *res1* to `F(z_1), F'(z_1)`
    by integrating the hypergeometric differential equation along a straight-line path.
    The evaluation points should be well-isolated from the singular points 0 and 1.

.. function:: void acb_hypgeom_2f1_series_direct(acb_poly_t res, const acb_poly_t a, const acb_poly_t b, const acb_poly_t c, const acb_poly_t z, int regularized, slong len, slong prec)

    Computes `F(z)` of the given power series truncated to length *len*, using
    direct summation of the hypergeometric series.

.. function:: void acb_hypgeom_2f1_direct(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec)

    Computes `F(z)` using direct summation of the hypergeometric series.

.. function:: void acb_hypgeom_2f1_transform(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, int which, slong prec)

.. function:: void acb_hypgeom_2f1_transform_limit(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, int which, slong prec)

    Computes `F(z)` using an argument transformation determined by the flag *which*.
    Legal values are 1 for `z/(z-1)`,
    2 for `1/z`, 3 for `1/(1-z)`, 4 for `1-z`, and 5 for `1-1/z`.

    The *limit* version assumes that *which* is not 1.
    If *which* is 2 or 3, it assumes that `b-a` represents an exact integer.
    If *which* is 4 or 5, it assumes that `c-a-b` represents an exact integer.
    In these cases, it computes the correct limit value.

.. function:: void acb_hypgeom_2f1_corner(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec)

    Computes `F(z)` near the corner cases `\exp(\pm \pi i \sqrt{3})`
    by analytic continuation.

.. function:: int acb_hypgeom_2f1_choose(const acb_t z)

    Chooses a method to compute the function based on the location of *z*
    in the complex plane. If the return value is 0, direct summation should be used.
    If the return value is 1 to 5, the transformation with this index in
    :func:`acb_hypgeom_2f1_transform` should be used.
    If the return value is 6, the corner case algorithm should be used.

.. function:: void acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec)

    Computes `F(z)` (or `\operatorname{\mathbf{F}}(z)` if *regularized* is set)
    using an automatic algorithm choice.

Orthogonal polynomials and functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_chebyshev_t(acb_t res, const acb_t n, const acb_t z, slong prec)

.. function:: void acb_hypgeom_chebyshev_u(acb_t res, const acb_t n, const acb_t z, slong prec)

    Computes the Chebyshev polynomial (or Chebyshev function) of first or second kind

    .. math ::

        T_n(z) = {}_2F_1\left(-n,n,\frac{1}{2},\frac{1-z}{2}\right)

    .. math ::

        U_n(z) = (n+1) {}_2F_1\left(-n,n+2,\frac{3}{2},\frac{1-z}{2}\right).

    The hypergeometric series definitions are only used for computation
    near the point 1. In general, trigonometric representations are used.
    For word-size integer *n*, :func:`acb_chebyshev_t_ui` and
    :func:`acb_chebyshev_u_ui` are called.

.. function:: void acb_hypgeom_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b, const acb_t z, slong prec)

    Computes the Jacobi polynomial (or Jacobi function)

    .. math ::

        P_n^{(a,b)}(z)=\frac{(a+1)_n}{\Gamma(n+1)} {}_2F_1\left(-n,n+a+b+1,a+1,\frac{1-z}{2}\right).

    For nonnegative integer *n*, this is a polynomial in *a*, *b* and *z*,
    even when the parameters are such that the hypergeometric series
    is undefined. In such cases, the polynomial is evaluated using
    direct methods.

.. function:: void acb_hypgeom_gegenbauer_c(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec)

    Computes the Gegenbauer polynomial (or Gegenbauer function)

    .. math ::

        C_n^{m}(z)=\frac{(2m)_n}{\Gamma(n+1)} {}_2F_1\left(-n,2m+n,m+\frac{1}{2},\frac{1-z}{2}\right).

    For nonnegative integer *n*, this is a polynomial in *m* and *z*,
    even when the parameters are such that the hypergeometric series
    is undefined. In such cases, the polynomial is evaluated using
    direct methods.

.. function:: void acb_hypgeom_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec)

    Computes the Laguerre polynomial (or Laguerre function)

    .. math ::

        L_n^{m}(z)=\frac{(m+1)_n}{\Gamma(n+1)} {}_1F_1\left(-n,m+1,z\right).

    For nonnegative integer *n*, this is a polynomial in *m* and *z*,
    even when the parameters are such that the hypergeometric series
    is undefined. In such cases, the polynomial is evaluated using
    direct methods.

    There are at least two incompatible ways to define the Laguerre function when
    *n* is a negative integer.  One possibility when `m = 0` is to define
    `L_{-n}^0(z) = e^z L_{n-1}^0(-z)`. Another possibility is to cover this
    case with the recurrence relation `L_{n-1}^m(z) + L_n^{m-1}(z) = L_n^m(z)`.
    Currently, we leave this case undefined (returning indeterminate).

.. function:: void acb_hypgeom_hermite_h(acb_t res, const acb_t n, const acb_t z, slong prec)

    Computes the Hermite polynomial (or Hermite function)

    .. math ::

        H_n(z) = 2^n \sqrt{\pi} \left(
            \frac{1}{\Gamma((1-n)/2)} {}_1F_1\left(-\frac{n}{2},\frac{1}{2},z^2\right)
            - 
            \frac{2z}{\Gamma(-n/2)} {}_1F_1\left(\frac{1-n}{2},\frac{3}{2},z^2\right)\right).

.. function:: void acb_hypgeom_legendre_p(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, slong prec)

    Sets *res* to the associated Legendre function of the first kind
    evaluated for degree *n*, order *m*, and argument *z*.
    When *m* is zero, this reduces to the Legendre polynomial `P_n(z)`.

    Many different branch cut conventions appear in the literature.
    If *type* is 0, the version

    .. math ::

        P_n^m(z) = \frac{(1+z)^{m/2}}{(1-z)^{m/2}}
            \mathbf{F}\left(-n, n+1, 1-m, \frac{1-z}{2}\right)

    is computed, and if *type* is 1, the alternative version

    .. math ::

        {\mathcal P}_n^m(z) = \frac{(z+1)^{m/2}}{(z-1)^{m/2}}
            \mathbf{F}\left(-n, n+1, 1-m, \frac{1-z}{2}\right).

    is computed. Type 0 and type 1 respectively correspond to
    type 2 and type 3 in *Mathematica* and *mpmath*.

.. function:: void acb_hypgeom_legendre_q(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, slong prec)

    Sets *res* to the associated Legendre function of the second kind
    evaluated for degree *n*, order *m*, and argument *z*.
    When *m* is zero, this reduces to the Legendre function `Q_n(z)`.

    Many different branch cut conventions appear in the literature.
    If *type* is 0, the version

    .. math ::

        Q_n^m(z) = \frac{\pi}{2 \sin(\pi m)}
            \left( \cos(\pi m) P_n^m(z) -
            \frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} P_n^{-m}(z)\right)

    is computed, and if *type* is 1, the alternative version

    .. math ::

        \mathcal{Q}_n^m(z) = \frac{\pi}{2 \sin(\pi m)} e^{\pi i m}
            \left( \mathcal{P}_n^m(z) -
            \frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} \mathcal{P}_n^{-m}(z)\right)

    is computed. Type 0 and type 1 respectively correspond to
    type 2 and type 3 in *Mathematica* and *mpmath*.

    When *m* is an integer, either expression is interpreted as a limit.
    We make use of the connection formulas [WQ3a]_, [WQ3b]_ and [WQ3c]_
    to allow computing the function even in the limiting case.
    (The formula [WQ3d]_ would be useful, but is incorrect in the lower
    half plane.)

    .. [WQ3a] http://functions.wolfram.com/07.11.26.0033.01
    .. [WQ3b] http://functions.wolfram.com/07.12.27.0014.01
    .. [WQ3c] http://functions.wolfram.com/07.12.26.0003.01
    .. [WQ3d] http://functions.wolfram.com/07.12.26.0088.01

.. function:: void acb_hypgeom_legendre_p_uiui_rec(acb_t res, ulong n, ulong m, const acb_t z, slong prec)

    For nonnegative integer *n* and *m*, uses recurrence relations to evaluate
    `(1-z^2)^{-m/2} P_n^m(z)` which is a polynomial in *z*.

.. function:: void acb_hypgeom_spherical_y(acb_t res, slong n, slong m, const acb_t theta, const acb_t phi, slong prec)

    Computes the spherical harmonic of degree *n*, order *m*,
    latitude angle *theta*, and longitude angle *phi*, normalized
    such that

    .. math ::

        Y_n^m(\theta, \phi) = \sqrt{\frac{2n+1}{4\pi} \frac{(n-m)!}{(n+m)!}} e^{im\phi} P_n^m(\cos(\theta)).

    The definition is extended to negative *m* and *n* by symmetry.
    This function is a polynomial in `\cos(\theta)` and `\sin(\theta)`.
    We evaluate it using :func:`acb_hypgeom_legendre_p_uiui_rec`.

