.. _acb-elliptic:

**acb_elliptic.h** -- elliptic integrals and functions of complex variables
===============================================================================


Carlson symmetric elliptic integrals
-------------------------------------------------------------------------------

Carlson symmetric forms are the preferred form of incomplete elliptic
integrals, due to their neat analytic structure and relatively
simple computation based on duplication theorems.
We largely follow the definitions and algorithms
in [Car1995]_ and chapter 19 in [NIST2012]_.

.. function:: void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)

    Evaluates the Carlson symmetric elliptic integral of the first kind

    .. math ::

        R_F(x,y,z) = \frac{1}{2}
            \int_0^{\infty} \frac{dt}{\sqrt{(t+x)(t+y)(t+z)}}

    where the square root extends continuously from the positive real axis.
    The function is well-defined for `x,y,z \notin (-\infty,0)`, and with
    at most one of `x,y,z` being zero.

    In general, one or more duplication steps are applied until
    `x,y,z` are close enough to use a multivariate Taylor polynomial
    of total degree 7.

    The special case `R_C(x, y) = R_F(x, y, y) = \frac{1}{2} \int_0^{\infty} (t+x)^{-1/2} (t+y)^{-1} dt`
    may be computed by
    setting *y* and *z* to the same variable.

    The *flags* parameter is reserved for future use and currently
    does nothing. Passing 0 results in default behavior.
