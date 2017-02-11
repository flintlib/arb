.. _acb-elliptic:

**acb_elliptic.h** -- elliptic integrals and functions of complex variables
===============================================================================

Warning: incomplete elliptic integrals have very complicated
branch structure when extended to complex variables.
For some functions in this module, branch cuts may be
artifacts of the evaluation algorithm rather than having
a natural mathematical justification.
The user should, accordingly, watch out for edge cases where the functions
implemented here may differ from other systems or literature.
There may also exist points where a function should be well-defined
but the implemented algorithm
fails to produce a finite result due to artificial internal singularities.

Carlson symmetric elliptic integrals
-------------------------------------------------------------------------------

Carlson symmetric forms are the preferred form of incomplete elliptic
integrals, due to their neat properties and relatively
simple computation based on duplication theorems.
We largely follow the definitions and algorithms
in [Car1995]_ and chapter 19 in [NIST2012]_.

.. function:: void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)

    Evaluates the Carlson symmetric elliptic integral of the first kind

    .. math ::

        R_F(x,y,z) = \frac{1}{2}
            \int_0^{\infty} \frac{dt}{\sqrt{(t+x)(t+y)(t+z)}}

    where the square root extends continuously from positive infinity.
    The integral is well-defined for `x,y,z \notin (-\infty,0)`, and with
    at most one of `x,y,z` being zero.
    When some parameters are negative real numbers, the function is
    still defined by analytic continuation.

    In general, one or more duplication steps are applied until
    `x,y,z` are close enough to use a multivariate Taylor polynomial
    of total degree 7.

    The special case `R_C(x, y) = R_F(x, y, y) = \frac{1}{2} \int_0^{\infty} (t+x)^{-1/2} (t+y)^{-1} dt`
    may be computed by
    setting *y* and *z* to the same variable.
    (This case is not yet handled specially, but might be optimized in
    the future.)

    The *flags* parameter is reserved for future use and currently
    does nothing. Passing 0 results in default behavior.

.. function:: void acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)

    Evaluates the Carlson symmetric elliptic integral of the third kind

    .. math ::

        R_J(x,y,z,p) = \frac{3}{2}
            \int_0^{\infty} \frac{dt}{(t+p)\sqrt{(t+x)(t+y)(t+z)}}.

    where the square root is taken continuously as in `R_J`.

    In general, one or more duplication steps are applied until
    `x,y,z,p` are close enough to use a multivariate Taylor polynomial
    of total degree 7.

    The duplication algorithm might not be correct for all possible
    combinations of complex variables, i.e. taking square roots
    during the computation might introduce spurious branch cuts.
    According to [Car1995]_, a sufficient (but not necessary) condition
    for correctness is that *x*, *y*, *z* have nonnegative
    real part and that *p* has positive real part.
    In other cases, the algorithm *may* still be correct, but the user
    should verify the results.

    The special case `R_D(x, y, z) = R_J(x, y, z, z)`
    may be computed by setting *z* and *p* to the same variable.
    This case is handled specially to avoid redundant arithmetic operations.
    In this case, the algorithm is correct for all *x*, *y* and *z*.

    The *flags* parameter is reserved for future use and currently
    does nothing. Passing 0 results in default behavior.

.. function:: void acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)

    Evaluates the Carlson symmetric elliptic integral of the second kind

    .. math ::

        R_G(x,y,z) = \frac{1}{4} \int_0^{\infty}
            \frac{t}{\sqrt{(t+x)(t+y)(t+z)}}
            \left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right) dt.

    The evaluation is done by expressing `R_G` in terms of `R_F` and `R_D`.
    There are no restrictions on the variables.

.. function:: void acb_elliptic_rc1(acb_t res, const acb_t x, slong prec)

    This helper function computes the special case
    `R_C(1, 1+x) = \operatorname{atan}(\sqrt{x})/\sqrt{x} = {}_2F_1(1,1/2,3/2,-x)`,
    which is needed in the evaluation of `R_J`.


