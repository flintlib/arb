.. _arb-hypgeom:

**arb_hypgeom.h** -- hypergeometric functions of real variables
==================================================================================

See :ref:`acb-hypgeom` for the general implementation of hypergeometric functions.

For convenience, this module provides versions of the same functions
for real variables
represented using :type:`arb_t` and :type:`arb_poly_t`.
Most methods are simple wrappers
around the complex versions,
but some of the functions in this module have been further optimized
specifically for real variables.

This module also provides certain functions exclusive to real variables,
such as functions for computing real roots of common special functions.

Generalized hypergeometric function
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_pfq(arb_t res, arb_srcptr a, slong p, arb_srcptr b, slong q, const arb_t z, int regularized, slong prec)

    Computes the generalized hypergeometric function `{}_pF_{q}(z)`,
    or the regularized version if *regularized* is set.

Confluent hypergeometric functions
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, slong prec)

    Computes the confluent hypergeometric limit function
    `{}_0F_1(a,z)`, or `\frac{1}{\Gamma(a)} {}_0F_1(a,z)` if *regularized*
    is set.

.. function:: void arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)

    Computes the confluent hypergeometric function
    `M(a,b,z) = {}_1F_1(a,b,z)`, or
    `\mathbf{M}(a,b,z) = \frac{1}{\Gamma(b)} {}_1F_1(a,b,z)` if *regularized* is set.

.. function:: void arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)

    Alias for :func:`arb_hypgeom_m`.

.. function:: void arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)

    Computes the confluent hypergeometric function `U(a,b,z)`.

Gauss hypergeometric function
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)

    Computes the Gauss hypergeometric function
    `{}_2F_1(a,b,c,z)`, or
    `\mathbf{F}(a,b,c,z) = \frac{1}{\Gamma(c)} {}_2F_1(a,b,c,z)`
    if *regularized* is set.

    Additional evaluation flags can be passed via the *regularized*
    argument; see :func:`acb_hypgeom_2f1` for documentation.

Error functions and Fresnel integrals
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_erf(arb_t res, const arb_t z, slong prec)

    Computes the error function `\operatorname{erf}(z)`.

.. function:: void _arb_hypgeom_erf_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_erf_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the error function of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_erfc(arb_t res, const arb_t z, slong prec)

    Computes the complementary error function
    `\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)`.
    This function avoids catastrophic cancellation for large positive *z*.

.. function:: void _arb_hypgeom_erfc_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_erfc_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the complementary error function of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_erfi(arb_t res, const arb_t z, slong prec)

    Computes the imaginary error function
    `\operatorname{erfi}(z) = -i\operatorname{erf}(iz)`.

.. function:: void _arb_hypgeom_erfi_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_erfi_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the imaginary error function of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, slong prec)

    Sets *res1* to the Fresnel sine integral `S(z)` and *res2* to
    the Fresnel cosine integral `C(z)`. Optionally, just a single function
    can be computed by passing *NULL* as the other output variable.
    The definition `S(z) = \int_0^z \sin(t^2) dt` is used if *normalized* is 0,
    and `S(z) = \int_0^z \sin(\tfrac{1}{2} \pi t^2) dt` is used if
    *normalized* is 1 (the latter is the Abramowitz & Stegun convention).
    `C(z)` is defined analogously.

.. function:: void _arb_hypgeom_fresnel_series(arb_ptr res1, arb_ptr res2, arb_srcptr z, slong zlen, int normalized, slong len, slong prec)

.. function:: void arb_hypgeom_fresnel_series(arb_poly_t res1, arb_poly_t res2, const arb_poly_t z, int normalized, slong len, slong prec)

    Sets *res1* to the Fresnel sine integral and *res2* to the Fresnel
    cosine integral of the power series *z*, truncated to length *len*.
    Optionally, just a single function can be computed by passing *NULL*
    as the other output variable.

Incomplete gamma and beta functions
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_gamma_upper(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)

    If *regularized* is 0, computes the upper incomplete gamma function
    `\Gamma(s,z)`.

    If *regularized* is 1, computes the regularized upper incomplete
    gamma function `Q(s,z) = \Gamma(s,z) / \Gamma(s)`.

    If *regularized* is 2, computes the generalized exponential integral
    `z^{-s} \Gamma(s,z) = E_{1-s}(z)` instead (this option is mainly
    intended for internal use; :func:`arb_hypgeom_expint` is the intended
    interface for computing the exponential integral).

.. function:: void _arb_hypgeom_gamma_upper_series(arb_ptr res, const arb_t s, arb_srcptr z, slong zlen, int regularized, slong n, slong prec)

.. function:: void arb_hypgeom_gamma_upper_series(arb_poly_t res, const arb_t s, const arb_poly_t z, int regularized, slong n, slong prec)

    Sets *res* to an upper incomplete gamma function where *s* is
    a constant and *z* is a power series, truncated to length *n*.
    The *regularized* argument has the same interpretation as in
    :func:`arb_hypgeom_gamma_upper`.

.. function:: void arb_hypgeom_gamma_lower(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)

    If *regularized* is 0, computes the lower incomplete gamma function
    `\gamma(s,z) = \frac{z^s}{s} {}_1F_1(s, s+1, -z)`.

    If *regularized* is 1, computes the regularized lower incomplete
    gamma function `P(s,z) = \gamma(s,z) / \Gamma(s)`.

    If *regularized* is 2, computes a further regularized lower incomplete
    gamma function `\gamma^{*}(s,z) = z^{-s} P(s,z)`.

.. function:: void _arb_hypgeom_gamma_lower_series(arb_ptr res, const arb_t s, arb_srcptr z, slong zlen, int regularized, slong n, slong prec)

.. function:: void arb_hypgeom_gamma_lower_series(arb_poly_t res, const arb_t s, const arb_poly_t z, int regularized, slong n, slong prec)

    Sets *res* to an lower incomplete gamma function where *s* is
    a constant and *z* is a power series, truncated to length *n*.
    The *regularized* argument has the same interpretation as in
    :func:`arb_hypgeom_gamma_lower`.

.. function:: void arb_hypgeom_beta_lower(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)

    Computes the (lower) incomplete beta function, defined by
    `B(a,b;z) = \int_0^z t^{a-1} (1-t)^{b-1}`,
    optionally the regularized incomplete beta function
    `I(a,b;z) = B(a,b;z) / B(a,b;1)`.

.. function:: void _arb_hypgeom_beta_lower_series(arb_ptr res, const arb_t a, const arb_t b, arb_srcptr z, slong zlen, int regularized, slong n, slong prec)

.. function:: void arb_hypgeom_beta_lower_series(arb_poly_t res, const arb_t a, const arb_t b, const arb_poly_t z, int regularized, slong n, slong prec)

    Sets *res* to the lower incomplete beta function `B(a,b;z)` (optionally
    the regularized version `I(a,b;z)`) where *a* and *b* are constants
    and *z* is a power series, truncating the result to length *n*.
    The underscore method requires positive lengths and does not support
    aliasing.

Exponential and trigonometric integrals
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_expint(arb_t res, const arb_t s, const arb_t z, slong prec)

    Computes the generalized exponential integral `E_s(z)`.

.. function:: void arb_hypgeom_ei(arb_t res, const arb_t z, slong prec)

    Computes the exponential integral `\operatorname{Ei}(z)`.

.. function:: void _arb_hypgeom_ei_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_ei_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the exponential integral of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_si(arb_t res, const arb_t z, slong prec)

    Computes the sine integral `\operatorname{Si}(z)`.

.. function:: void _arb_hypgeom_si_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_si_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the sine integral of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_ci(arb_t res, const arb_t z, slong prec)

    Computes the cosine integral `\operatorname{Ci}(z)`.
    The result is indeterminate if `z < 0` since the value of the function would be complex.

.. function:: void _arb_hypgeom_ci_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_ci_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the cosine integral of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_shi(arb_t res, const arb_t z, slong prec)

    Computes the hyperbolic sine integral `\operatorname{Shi}(z) = -i \operatorname{Si}(iz)`.

.. function:: void _arb_hypgeom_shi_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_shi_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the hyperbolic sine integral of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_chi(arb_t res, const arb_t z, slong prec)

    Computes the hyperbolic cosine integral `\operatorname{Chi}(z)`.
    The result is indeterminate if `z < 0` since the value of the function would be complex.

.. function:: void _arb_hypgeom_chi_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_chi_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

    Computes the hyperbolic cosine integral of the power series *z*,
    truncated to length *len*.

.. function:: void arb_hypgeom_li(arb_t res, const arb_t z, int offset, slong prec)

    If *offset* is zero, computes the logarithmic integral
    `\operatorname{li}(z) = \operatorname{Ei}(\log(z))`.

    If *offset* is nonzero, computes the offset logarithmic integral
    `\operatorname{Li}(z) = \operatorname{li}(z) - \operatorname{li}(2)`.

    The result is indeterminate if `z < 0` since the value of the function would be complex.

.. function:: void _arb_hypgeom_li_series(arb_ptr res, arb_srcptr z, slong zlen, int offset, slong len, slong prec)

.. function:: void arb_hypgeom_li_series(arb_poly_t res, const arb_poly_t z, int offset, slong len, slong prec)

    Computes the logarithmic integral (optionally the offset version)
    of the power series *z*, truncated to length *len*.

Bessel functions
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_bessel_j(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes the Bessel function of the first kind `J_{\nu}(z)`.

.. function:: void arb_hypgeom_bessel_y(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes the Bessel function of the second kind `Y_{\nu}(z)`.

.. function:: void arb_hypgeom_bessel_jy(arb_t res1, arb_t res2, const arb_t nu, const arb_t z, slong prec)

    Sets *res1* to `J_{\nu}(z)` and *res2* to `Y_{\nu}(z)`, computed
    simultaneously.

.. function:: void arb_hypgeom_bessel_i(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes the modified Bessel function of the first kind
    `I_{\nu}(z) = z^{\nu} (iz)^{-\nu} J_{\nu}(iz)`.

.. function:: void arb_hypgeom_bessel_i_scaled(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes the function `e^{-z} I_{\nu}(z)`.

.. function:: void arb_hypgeom_bessel_k(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes the modified Bessel function of the second kind `K_{\nu}(z)`.

.. function:: void arb_hypgeom_bessel_k_scaled(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes the function `e^{z} K_{\nu}(z)`.

Airy functions
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_airy(arb_t ai, arb_t ai_prime, arb_t bi, arb_t bi_prime, const arb_t z, slong prec)

    Computes the Airy functions `(\operatorname{Ai}(z), \operatorname{Ai}'(z), \operatorname{Bi}(z), \operatorname{Bi}'(z))`
    simultaneously. Any of the four function values can be omitted by passing
    *NULL* for the unwanted output variables, speeding up the evaluation.

.. function:: void arb_hypgeom_airy_jet(arb_ptr ai, arb_ptr bi, const arb_t z, slong len, slong prec)

    Writes to *ai* and *bi* the respective Taylor expansions of the Airy functions
    at the point *z*, truncated to length *len*.
    Either of the outputs can be *NULL* to avoid computing that function.
    The variable *z* is not allowed to be aliased with the outputs.
    To simplify the implementation, this method does not compute the
    series expansions of the primed versions directly; these are
    easily obtained by computing one extra coefficient and differentiating
    the output with :func:`_arb_poly_derivative`.

.. function:: void _arb_hypgeom_airy_series(arb_ptr ai, arb_ptr ai_prime, arb_ptr bi, arb_ptr bi_prime, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_airy_series(arb_poly_t ai, arb_poly_t ai_prime, arb_poly_t bi, arb_poly_t bi_prime, const arb_poly_t z, slong len, slong prec)

    Computes the Airy functions evaluated at the power series *z*,
    truncated to length *len*. As with the other Airy methods, any of the
    outputs can be *NULL*.

.. function:: void arb_hypgeom_airy_zero(arb_t a, arb_t a_prime, arb_t b, arb_t b_prime, const fmpz_t n, slong prec)

    Computes the *n*-th real zero `a_n`, `a'_n`, `b_n`, or `b'_n`
    for the respective Airy function or Airy function derivative.
    Any combination of the four output variables can be *NULL*.
    The zeros are indexed by increasing magnitude, starting with
    `n = 1` to follow the convention in the literature.
    An index *n* that is not positive is invalid input.
    The implementation uses asymptotic expansions for the zeros
    [PS1991]_ together with the interval Newton method for refinement.

Coulomb wave functions
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_coulomb(arb_t F, arb_t G, const arb_t l, const arb_t eta, const arb_t z, slong prec)

    Writes to *F*, *G* the values of the respective
    Coulomb wave functions `F_{\ell}(\eta,z)` and `G_{\ell}(\eta,z)`.
    Either of the outputs can be *NULL*.

.. function:: void arb_hypgeom_coulomb_jet(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, const arb_t z, slong len, slong prec)

    Writes to *F*, *G* the respective Taylor expansions of the
    Coulomb wave functions at the point *z*, truncated to length *len*.
    Either of the outputs can be *NULL*.

.. function:: void _arb_hypgeom_coulomb_series(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, arb_srcptr z, slong zlen, slong len, slong prec)

.. function:: void arb_hypgeom_coulomb_series(arb_poly_t F, arb_poly_t G, const arb_t l, const arb_t eta, const arb_poly_t z, slong len, slong prec)

    Computes the Coulomb wave functions evaluated at the power series *z*,
    truncated to length *len*. Either of the outputs can be *NULL*.

Orthogonal polynomials and functions
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_chebyshev_t(arb_t res, const arb_t nu, const arb_t z, slong prec)

.. function:: void arb_hypgeom_chebyshev_u(arb_t res, const arb_t nu, const arb_t z, slong prec)

.. function:: void arb_hypgeom_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t z, slong prec)

.. function:: void arb_hypgeom_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t z, slong prec)

.. function:: void arb_hypgeom_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t z, slong prec)

.. function:: void arb_hypgeom_hermite_h(arb_t res, const arb_t nu, const arb_t z, slong prec)

    Computes Chebyshev, Jacobi, Gegenbauer, Laguerre or Hermite polynomials,
    or their extensions to non-integer orders.

.. function:: void arb_hypgeom_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, slong prec)

.. function:: void arb_hypgeom_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, slong prec)

    Computes Legendre functions of the first and second kind.
    See :func:`acb_hypgeom_legendre_p` and :func:`acb_hypgeom_legendre_q`
    for definitions.

.. function:: void arb_hypgeom_legendre_p_ui_deriv_bound(mag_t dp, mag_t dp2, ulong n, const arb_t x, const arb_t x2sub1)

    Sets *dp* to an upper bound for `P'_n(x)` and *dp2* to an upper
    bound for `P''_n(x)` given *x* assumed to represent a real
    number with `|x| \le 1`. The variable *x2sub1* must contain
    the precomputed value `1-x^2` (or `x^2-1`). This method is used
    internally to bound the propagated error for Legendre polynomials.

.. function:: void arb_hypgeom_legendre_p_ui_zero(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong K, slong prec)

.. function:: void arb_hypgeom_legendre_p_ui_one(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong K, slong prec)

.. function:: void arb_hypgeom_legendre_p_ui_asymp(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong K, slong prec)

.. function:: void arb_hypgeom_legendre_p_rec(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)

.. function:: void arb_hypgeom_legendre_p_ui(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)

    Evaluates the ordinary Legendre polynomial `P_n(x)`. If *res_prime* is
    non-NULL, simultaneously evaluates the derivative `P'_n(x)`.

    The overall algorithm is described in [JM2018]_.

    The versions *zero*, *one* respectively use the hypergeometric series
    expansions at `x = 0` and `x = 1` while the *asymp* version uses an
    asymptotic series on `(-1,1)` intended for large *n*. The parameter *K*
    specifies the exact number of expansion terms to use (if the series
    expansion truncated at this point does not give the exact polynomial,
    an error bound is computed automatically).
    The asymptotic expansion with error bounds is given in [Bog2012]_.
    The *rec* version uses the forward recurrence implemented using
    fixed-point arithmetic; it is only intended for the interval `(-1,1)`,
    moderate *n* and modest precision.

    The default version attempts to choose the best algorithm automatically.
    It also estimates the amount of cancellation in the hypergeometric series
    and increases the working precision to compensate, bounding the
    propagated error using derivative bounds.

.. function:: void arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, ulong n, ulong k, slong prec)

    Sets *res* to the *k*-th root of the Legendre polynomial `P_n(x)`.
    We index the roots in decreasing order

    .. math ::

        1 > x_0 > x_1 > \ldots > x_{n-1} > -1

    (which corresponds to ordering the roots of `P_n(\cos(\theta))`
    in order of increasing `\theta`).
    If *weight* is non-NULL, it is set to the weight corresponding
    to the node `x_k` for Gaussian quadrature on `[-1,1]`.
    Note that only `\lceil n / 2 \rceil` roots need to be computed,
    since the remaining roots are given by `x_k = -x_{n-1-k}`.

    We compute an enclosing interval using an asymptotic approximation followed
    by some number of Newton iterations, using the error bounds given
    in [Pet1999]_. If very high precision is requested, the root is
    subsequently refined using interval Newton steps with doubling working
    precision.

Dilogarithm
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_dilog(arb_t res, const arb_t z, slong prec)

    Computes the dilogarithm `\operatorname{Li}_2(z)`.

Hypergeometric sequences
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_central_bin_ui(arb_t res, ulong n, slong prec)

    Computes the central binomial coefficient `{2n \choose n}`.

