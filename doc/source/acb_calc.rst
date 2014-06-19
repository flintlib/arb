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

        int func(acb_ptr out, const acb_t inp, void * param, long order, long prec)

    implementing a univariate complex function `f(x)`.
    When called, *func* should write to *out* the first *order*
    coefficients in the Taylor series expansion of `f(x)` at the point *inp*,
    evaluated at a precision of *prec* bits.
    The *param* argument may be used to pass through
    additional parameters to the function.
    The return value is reserved for future use as an
    error code. It can be assumed that *out* and *inp* are not
    aliased and that *order* is positive.


Bounds
-------------------------------------------------------------------------------

.. function:: void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param, const acb_t x, const arb_t radius, long maxdepth, long prec)

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

Integration
-------------------------------------------------------------------------------

.. function:: int acb_calc_integrate_taylor(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const arf_t inner_radius, const arf_t outer_radius, long accuracy_goal, long prec)

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

