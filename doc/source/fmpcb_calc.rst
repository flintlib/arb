.. _fmpcb-calc:

**fmpcb_calc.h** -- calculus with complex-valued functions
===============================================================================

This module provides functions for operations of calculus
over the complex numbers (intended to include root-finding,
integration, and so on).



Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpcb_calc_func_t

    Typedef for a pointer to a function with signature

    .. code ::

        int func(fmpcb_ptr out, const fmpcb_t inp, void * param, long order, long prec)

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

.. function :: void fmpcb_calc_cauchy_bound(fmprb_t bound, fmpcb_calc_func_t func, void * param, const fmpcb_t x, const fmprb_t radius, long maxdepth, long prec)

    Sets *bound* to a ball containing the value of the integral

    .. math ::

        b = \frac{1}{2 \pi r} \oint_{|t-x| = r} |f(t)| dt

    where *f* is specified by (*func*, *param*) and *r* is given by *radius*.
    The integral is computed using a simple step sum.
    The integration range is subdivided until the order of magnitude of *b*
    can be determined (i.e. its error bound is smaller than its midpoint),
    or until the step length has been cut in half *maxdepth* times.
    This function is currently implemented completely naively, and
    repeatedly subdivides the whole integration range instead of
    performing adaptive subdivisions.

