.. _acb-calc:

**acb_calc.h** -- calculus with complex-valued functions
===============================================================================

This module provides functions for operations of calculus
over the complex numbers (intended to include root-finding,
integration, and so on).

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: acb_calc_func_t

    Typedef for a pointer to a function with signature::

        int func(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)

    implementing a univariate complex function `f(z)`.
    The *param* argument may be used to pass through
    additional parameters to the function.
    The return value is reserved for future use as an
    error code. It can be assumed that *out* and *inp* are not aliased.

    When called with *order* = 0, *func* should write to *out* the value
    of `f(z)` at the point *inp*, evaluated at a precision of *prec* bits.
    In this case, *f* can be an arbitrary complex function, which
    may have branch cuts or even be non-holomorphic.

    When called with *order* = *n* for `n \ge 1`, *func* should write to
    *out* the first *n* coefficients in the Taylor series expansion of `f(z)`
    at the point *inp*, evaluated at a precision of *prec* bits.
    In this case, the implementation of *func* must verify that *f*
    is holomorphic on the complex interval defined by *z*, and set the
    coefficients in *out* to non-finite values otherwise.

    For algorithms that do not rely on derivatives, *func* will always
    get called with *order* = 0 or *order* = 1, in which case the user
    only needs to implement evaluation of the direct function value `f(z)`
    (without derivatives). With *order* = 1, *func* must verify
    holomorphicity (unlike the *order* = 0 case).

    If *f* is built from field operations and meromorphic functions, then
    no special action is necessary when *order* is positive
    since division by zero or evaluation
    of builtin functions at poles automatically produces infinite enclosures.
    However, manual action is needed for bounded functions with branch cuts.
    For example, when evaluating `\sqrt{z}`, the output must be set to
    an non-finite value if `z` overlaps with the branch cut `[-\infty,0]`.
    The easiest solution is to use holomorphy-testing versions of atomic
    functions.

    For example, here we define a piecewise holomorphic extension
    of the function
    `f(z) = \sqrt{\lfloor z \rfloor}` (for simplicity, without implementing
    derivatives)::

        /* Floor function on R extended to a piecewise holomorphic function in
           vertical strips. */
        void holomorphic_floor(acb_t res, const acb_t z, int holomorphic, slong prec)
        {
            if (!acb_is_finite(z) || (holomorphic && arb_contains_int(acb_realref(z))))
            {
                acb_indeterminate(res);
            }
            else
            {
                arb_floor(acb_realref(res), acb_realref(z), prec);
                arb_set_round(acb_imagref(res), acb_imagref(z), prec);
            }
        }

        /* Square root function on C with detection of the branch cut. */
        void holomorphic_sqrt(acb_t res, const acb_t z, int holomorphic, slong prec)
        {
            if (!acb_is_finite(z) || (holomorphic &&
                arb_contains_zero(acb_imagref(z)) &&
                arb_contains_nonpositive(acb_realref(z))))
            {
                acb_indeterminate(res);
            }
            else
            {
                acb_sqrt(res, z, prec);
            }
        }

        int func(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)
        {
            if (order > 1) flint_abort();  /* derivatives not implemented */

            holomorphic_floor(out, inp, order != 0, prec);
            holomorphic_sqrt(out, out, order != 0, prec);
            return 0;
        }


Integration
-------------------------------------------------------------------------------

.. function:: void acb_calc_integrate(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, slong goal, const mag_t tol, slong deg_limit, slong eval_limit, slong depth_limit, int flags, slong prec)

    Computes a rigorous enclosure of the integral

    .. math ::

        I = \int_a^b f(t) dt

    where *f* is specified by (*func*, *param*), following a straight-line
    path between the complex numbers *a* and *b*.
    For finite results, *a*, *b* must be finite and *f* must be bounded
    on the path of integration.
    To compute improper integrals, the user should therefore truncate the path
    of integration manually (or make a regularizing change of variables,
    if possible).

    By default, the integrand *func* will only be called with *order* = 0
    or *order* = 1; that is, derivatives are not required.

    - The integrand will be called with *order* = 0 to evaluate *f*
      normally on the integration path (either at a single point
      or on a subinterval). In this case, *f* is treated as a pointwise defined
      function and can have arbitrary discontinuities.

    - The integrand will be called with *order* = 1 to evaluate *f*
      on a domain surrounding a segment of the integration path for the purpose
      of bounding the error of a quadrature formula. In this case, *func* must
      verify that *f* is holomorphic on this domain (and output a non-finite
      value if it is not).

    The integration algorithm combines direct interval enclosures,
    Gauss-Legendre quadrature where *f* is holomorphic,
    and adaptive subdivision. This strategy supports integrands with
    discontinuities while providing exponential convergence for typical
    piecewise holomorphic integrands.

    On each subinterval, the algorithm attempts to
    achieve a relative error of `2^{-\text{goal}}` or an absolute error of
    *tol*, whichever is larger. The parameters *goal* and *tol* are only
    guidelines; the cumulative error may be larger than both the prescribed
    absolute and relative error goals, depending on the number of
    subdivisions, cancellation between segments of the integral, and numerical
    errors in the evaluation of the integrand.

    The following parameters control the integration.

    - *prec* - working precision. This is the working precision used to
      evaluate the integrand and manipulate interval endpoints.

    - *goal* - relative accuracy goal as a nonnegative number of bits, i.e.
      target a relative error of `2^{-\text{goal}}`. This parameter can
      simply be set to *prec* (or a slightly smaller value) in most situations.

    - *tol* - absolute tolerance goal (specified as a :type:`mag_t`).
      In general, *tol* should be set to about `2^{-\text{goal}}` times
      an estimate for the magnitude of the integral. For typical
      integrals which will have magnitude within a few orders of magnitude
      of unity, `\text{tol} = 2^{-\text{goal}}` works well enough.
      Of course, if the integral has very small magnitude, the absolute
      tolerance must be set to a small value to get high relative accuracy.

      If the integral has very large magnitude, setting the absolute
      tolerance to a corresponding large value is recommended for best
      performance, but it is not necessary for convergence since the absolute
      tolerance is increased automatically during the execution of the
      algorithm if the partial integrals are found to have larger error.

      Setting *tol* to 0 is allowed and forces use of a relative instead of an
      absolute tolerance goal, which can be handy for exponentially small or
      large functions of unknown magnitude. It is recommended to avoid this
      solution if possible since the algorithm might need many extra
      subdivisions to determine an initial scale; if the approximate
      magnitude can be estimated by some external means (for example if
      a midpoint-width or endpoint-width estimate is known to be accurate),
      the caller should set *tol* to `2^{-\text{goal}}`
      times such an external estimate for best performance.

    - *deg_limit* - maximum quadrature degree for each subinterval.
      If a zero or negative value is provided, the limit is set to a default
      value which currently equals `0.5 \cdot \text{goal} + 10` for
      Gauss-Legendre quadrature.

      A higher quadrature degree can be beneficial for functions that
      are holomorphic on a large domain around the integration path
      and yet behave irregularly, such as entire functions with a high
      amount of oscillation. The drawback of increasing the degree is that
      the precomputation time for quadrature nodes increases.

    - *eval_limit* - maximum number of function evaluations.
      If a zero or negative value is provided, the limit is set to a default
      value which currently equals `1000 \cdot \text{prec} + \text{prec}^2`.

      This is the main parameter used to limit the amount of work before
      aborting due to possible slow convergence or non-convergence.
      (This limit is only taken as a rough guideline, and the actual number of
      function evaluations may be slightly higher depending on the
      actual subdivisions.)
      A lower limit allows aborting faster. A higher limit may be needed
      for integrands with many discontinuities or many singularities
      close to the integration path.

    - *depth_limit* - maximum subdivision depth.
      If a zero or negative value is provided, the limit is set to the
      default value `2 \cdot \text{prec}`.
      This limits the number of recursive bisections around any single point.
      There is typically no reason to use a non-default value.

    - *flags* - additional options

        *ACB_CALC_VERBOSE*          - print some information

        *ACB_CALC_VERY_VERBOSE*      - print even more information

.. function:: slong acb_calc_integrate_gl_auto_deg(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const mag_t tol, slong deg_limit, int flags, slong prec)

    Attempts to compute `I = \int_a^b f(t) dt` using a single application
    of Gauss-Legendre quadrature with automatic determination of the
    quadrature degree so that the error is smaller than *tol*.
    Returns a positive integer *n* indicating that the integral has
    been evaluated successfully, or returns 0 if the tolerance could not be met.

    For the interval `[-1,1]`, the error of the *n*-point Gauss-Legendre
    rule is bounded by

    .. math ::

        \left| I - \sum_{k=0}^{n-1} w_k f(x_k) \right| \le \frac{64 M}{15 (\rho-1) \rho^{2n-1}}

    if `f` is holomorphic with `|f(z)| \le M` inside the ellipse *E*
    with foci `\pm 1` and semiaxes
    `X` and `Y = \sqrt{X^2 - 1}` such that `\rho = X + Y`
    with `\rho > 1`
    (See Trefethen, "Is Gauss Quadrature Better than Clenshaw-Curtis?").

    For an arbitrary interval, we use `\int_a^b f(t) dt = \int_{-1}^1 g(t) dt`
    where `g(t) = \Delta f(\Delta t + m)`,
    `\Delta = \tfrac{1}{2}(b-a)`, `m = \tfrac{1}{2}(a+b)`.
    With `I = [\pm X] + [\pm Y]i`, this means that we evaluate
    `\Delta f(\Delta I + m)` to get the bound `M`.
    (An improvement would be to reduce the wrapping effect of rotating the
    ellipse when the path is not rectilinear).

    We search for an `X` that makes the error small by trying steps `2^{2^k}`.
    Larger `X` will give smaller `1 / \rho^{2n-1}` but larger `M`. If we try
    successive larger values of `k`, we can abort when `M = \infty`
    since this either means that we have hit a singularity or a branch cut or
    that overestimation in the evaluation of `f` is becoming too severe.

Integration (old)
-------------------------------------------------------------------------------

.. function:: void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param, const acb_t x, const arb_t radius, slong maxdepth, slong prec)

    Sets *bound* to a ball containing the value of the integral

    .. math ::

        C(x,r) = \frac{1}{2 \pi r} \oint_{|z-x| = r} |f(z)| dz
               = \int_0^1 |f(x+re^{2\pi i t})| dt

    where *f* is specified by (*func*, *param*) and *r* is given by *radius*.
    The integral is computed using a simple step sum.
    The integration range is subdivided until the order of magnitude of *b*
    can be determined (i.e. its error bound is smaller than its midpoint),
    or until the step length has been cut in half *maxdepth* times.
    This function is currently implemented completely naively, and
    repeatedly subdivides the whole integration range instead of
    performing adaptive subdivisions.

.. function:: int acb_calc_integrate_taylor(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const arf_t inner_radius, const arf_t outer_radius, slong accuracy_goal, slong prec)

    Computes the integral

    .. math ::

        I = \int_a^b f(t) dt

    where *f* is specified by (*func*, *param*), following a straight-line
    path between the complex numbers *a* and *b* which both must be finite.

    The integral is approximated by piecewise centered Taylor polynomials.
    Rigorous truncation error bounds are calculated using the Cauchy integral
    formula. More precisely, if the Taylor series of *f* centered at the point
    *m* is `f(m+x) = \sum_{n=0}^{\infty} a_n x^n`, then

    .. math ::

        \int f(m+x) = \left( \sum_{n=0}^{N-1} a_n \frac{x^{n+1}}{n+1} \right)
                  + \left( \sum_{n=N}^{\infty} a_n \frac{x^{n+1}}{n+1} \right).

    For sufficiently small *x*, the second series converges and its
    absolute value is bounded by

    .. math ::

        \sum_{n=N}^{\infty} \frac{C(m,R)}{R^n} \frac{|x|^{n+1}}{N+1}
            = \frac{C(m,R) R x}{(R-x)(N+1)} \left( \frac{x}{R} \right)^N.

    It is required that any singularities of *f* are
    isolated from the path of integration by a distance strictly
    greater than the positive value *outer_radius* (which is the integration
    radius used for the Cauchy bound). Taylor series step lengths are
    chosen so as not to
    exceed *inner_radius*, which must be strictly smaller than *outer_radius*
    for convergence. A smaller *inner_radius* gives more rapid convergence
    of each Taylor series but means that more series might have to be used.
    A reasonable choice might be to set *inner_radius* to half the value of
    *outer_radius*, giving roughly one accurate bit per term.

    The truncation point of each Taylor series is chosen so that the absolute
    truncation error is roughly `2^{-p}` where *p* is given by *accuracy_goal*
    (in the future, this might change to a relative accuracy).
    Arithmetic operations and function
    evaluations are performed at a precision of *prec* bits. Note that due
    to accumulation of numerical errors, both values may have to be set
    higher (and the endpoints may have to be computed more accurately)
    to achieve a desired accuracy.

    This function chooses the evaluation points uniformly rather
    than implementing adaptive subdivision.

