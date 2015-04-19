.. _acb-hypgeom:

**acb_hypgeom.h** -- hypergeometric functions in the complex numbers
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

    {}_pH_{q}(a_0,\ldots,a_{p-1}; b_0 \ldots b_{q-1}; z) = {}_{p+1}F_{q}(a_0,\ldots,a_{p-1},1; b_0 \ldots b_{q-1}; z) = \sum_{k=0}^{\infty} T(k)

For the conventional generalized hypergeometric function
`{}_pF_{q}`, compute  `{}_pH_{q+1}` with the explicit parameter `b_q = 1`,
or remove a 1 from the `a_i` parameters if there is one.

.. function:: void acb_hypgeom_pfq_bound_factor(mag_t C, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, ulong n)

    Computes a factor *C* such that

    .. math ::

        \left|\sum_{k=n}^{\infty} T(k)\right| \le C |T(n)|.

    We check that `\operatorname{Re}(b+n) > 0` for all lower
    parameters *b*. If this does not hold, *C* is set to infinity.
    Otherwise, we cancel out pairs of parameters
    `a` and `b` against each other. We have

    .. math ::

        \left|\frac{a+k}{b+k}\right| = \left|1 + \frac{a-b}{b+k}\right| \le 1 + \frac{|a-b|}{|b+n|}

    and

    .. math ::

        \left|\frac{1}{b+k}\right| \le \frac{1}{|b+n|}

    for all `k \ge n`. This gives us a constant *D* such that
    `T(k+1) \le D T(k)` for all `k \ge n`.
    If `D \ge 1`, we set *C* to infinity. Otherwise, we take
    `C = \sum_{k=0}^{\infty} D^k = (1-D)^{-1}`.

    As currently implemented, the bound becomes infinite when `n` is
    too small, even if the series converges.

.. function:: long acb_hypgeom_pfq_choose_n(acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long prec)

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

.. function:: void acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)

.. function:: void acb_hypgeom_pfq_sum_rs(acb_t s, acb_t t, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)

.. function:: void acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)

    Computes `s = \sum_{k=0}^{n-1} T(k)` and `t = T(n)`.
    Does not allow aliasing between input and output variables.
    We require `n \ge 0`.

    The *forward* version computes the sum using forward
    recurrence.

    The *rs* version computes the sum in reverse order
    using rectangular splitting. It only computes a
    magnitude bound for the value of *t*.

    The default version automatically chooses an algorithm
    depending on the inputs.

.. function:: void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)

    Computes

    .. math ::

        {}_pH_{q}(z)
            = \sum_{k=0}^{\infty} T(k)
            = \sum_{k=0}^{n-1} T(k) + \varepsilon

    directly from the defining series, including a rigorous bound for
    the truncation error `\varepsilon` in the output.

    If  `n < 0`, this function chooses a number of terms automatically
    using :func:`acb_hypgeom_pfq_choose_n`.

.. function:: void acb_hypgeom_pfq_series_direct(acb_poly_t res, const acb_poly_struct * a, long p, const acb_poly_struct * b, long q, const acb_poly_t z, int regularized, long n, long len, long prec)

    Computes `{}_pH_{q}(z)` directly using the defining series, given
    parameters and argument that are power series.
    The result is a power series of length *len*.

    An error bound is computed automatically as a function of the number
    of terms *n*. If `n < 0`, the number of terms is chosen
    automatically.

    If *regularized* is set, the regularized hypergeometric function
    is computed instead.

Asymptotic series
-------------------------------------------------------------------------------

Let `U(a,b,z)` denote the confluent hypergeometric function of the second
kind with the principal branch cut, and
let `U^{*} = z^a U(a,b,z)`.
For all `z \ne 0` and `b \notin \mathbb{Z}` (but valid for all `b` as a limit),
we have (DLMF 13.2.42)

.. math ::

    U(a,b,z)
        = \frac{\Gamma(1-b)}{\Gamma(a-b+1)} M(a,b,z)
        + \frac{\Gamma(b-1)}{\Gamma(a)} z^{1-b} M(a-b+1,2-b,z).

Moreover, for all `z \ne 0` we have

.. math ::

    \frac{{}_1F_1(a,b,z)}{\Gamma(b)}
        = \frac{(-z)^{-a}}{\Gamma(b-a)} U^{*}(a,b,z)
        + \frac{z^{a-b} e^z}{\Gamma(a)} U^{*}(b-a,b,-z)

which is equivalent to DLMF 13.2.41 (but simpler in form).

We have the asymptotic expansion

.. math ::

    U^{*}(a,b,z) \sim {}_2F_0(a, a-b+1, -1/z)

where `{}_2F_0(a,b,z)` denotes a formal hypergeometric series, i.e.

.. math ::

    U^{*}(a,b,z) = \sum_{k=0}^{n-1} \frac{(a)_k (a-b+1)_k}{k! (-z)^k} + \varepsilon_n(z).

The error term `\varepsilon_n(z)` is bounded according to DLMF 13.7.
A case distinction is made depending on whether `z` lies in one
of three regions which we index by `R`.
Our formula for the error bound increases with the value of `R`, so we
can always choose the larger out of two indices if `z` lies in
the union of two regions.

Let `r = |b-2a|`.
If `\operatorname{Re}(z) \ge r`, set `R = 1`.
Otherwise, if `\operatorname{Im}(z) \ge r` or `\operatorname{Re}(z) \ge 0 \land |z| \ge r`, set `R = 2`.
Otherwise, if `|z| \ge 2r`, set `R = 3`.
Otherwise, the bound is infinite.
If the bound is finite, we have

.. math ::

    |\varepsilon_n(z)| \le 2 \alpha C_n \left|\frac{(a)_n (a-b+1)_n}{n! z^n} \right| \exp(2 \alpha \rho C_1 / |z|)

in terms of the following auxiliary quantities

.. math ::

    \sigma = |(b-2a)/z|

    C_n = \begin{cases}
    1                              & \text{if } R = 1 \\
    \chi(n)                        & \text{if } R = 2 \\
    (\chi(n) + \rho \nu^2 n) \nu^n & \text{if } R = 3
    \end{cases}

    \nu = \left(\tfrac{1}{2} + \tfrac{1}{2}\sqrt{1-4\sigma^2}\right)^{-1/2} \le 1 + 2 \sigma^2

    \chi(n) = \sqrt{\pi} \Gamma(\tfrac{1}{2}n+1) / \Gamma(\tfrac{1}{2} n + \tfrac{1}{2})

    \sigma' = \begin{cases}
    \sigma & \text{if } R \ne 3 \\
    \nu \sigma & \text{if } R = 3
    \end{cases}

    \alpha = (1 - \sigma')^{-1}

    \rho = \tfrac{1}{2} |2a^2-2ab+b| + \sigma' (1+ \tfrac{1}{4} \sigma') (1-\sigma')^{-2}

.. function:: void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, long n, long prec)

    Sets *res* to `U^{*}(a,b,z)` computed using *n* terms of the asymptotic series,
    with a rigorous bound for the error included in the output.
    We require `n \ge 0`.

.. function:: int acb_hypgeom_u_use_asymp(const acb_t z, long prec)

    Heuristically determines whether the asymptotic series can be used
    to evaluate `U(a,b,z)` to *prec* accurate bits (assuming that *a* and *b*
    are small).

Confluent hypergeometric functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_u_1f1_series(acb_poly_t res, const acb_poly_t a, const acb_poly_t b, const acb_poly_t z, long len, long prec)

    Computes `U(a,b,z)` as a power series truncated to length *len*,
    given `a, b, z \in \mathbb{C}[[x]]`.
    If `b[0] \in \mathbb{Z}`, it computes one extra derivative and removes
    the singularity (it is then assumed that `b[1] \ne 0`).
    As currently implemented, the output is indeterminate if `b` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)

    Computes `U(a,b,z)` as a sum of two convergent hypergeometric series.
    If `b \in \mathbb{Z}`, it computes
    the limit value via :func:`acb_hypgeom_u_1f1_series`.
    As currently implemented, the output is indeterminate if `b` is nonexact
    and contains an integer.

.. function:: void acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)

    Computes `U(a,b,z)` using an automatic algorithm choice. The
    function :func:`acb_hypgeom_u_asymp` is used
    if `a` or `a-b+1` is a nonpositive integer (in which
    case the asymptotic series terminates), or if *z* is sufficiently large.
    Otherwise :func:`acb_hypgeom_u_1f1` is used.

.. function:: void acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)

.. function:: void acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)

.. function:: void acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)

    Computes the confluent hypergeometric function
    `M(a,b,z) = {}_1F_1(a,b,z)`, or
    `\mathbf{M}(a,b,z) = \frac{1}{\Gamma(b)} {}_1F_1(a,b,z)` if *regularized*
    is set.

The error function
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_erf_asymp(acb_t res, const acb_t z, long prec, long prec2)

.. function:: void acb_hypgeom_erf(acb_t res, const acb_t z, long prec)

    Computes the error function respectively using

    .. math ::

        \operatorname{erf}(z) = \frac{2z}{\sqrt{\pi}}
            {}_1F_1(\tfrac{1}{2}, \tfrac{3}{2}, -z^2)

        \operatorname{erf}(z) = \frac{2z e^{-z^2}}{\sqrt{\pi}}
            {}_1F_1(1, \tfrac{3}{2}, z^2)

        \operatorname{erf}(z) = \frac{z}{\sqrt{z^2}}
            \left(1 - \frac{e^{-z^2}}{\sqrt{\pi}}
            U(\tfrac{1}{2}, \tfrac{1}{2}, z^2)\right).

    and an automatic algorithm choice. The *asymp* version takes a second
    precision to use for the *U* term.

.. function:: void acb_hypgeom_erfc(acb_t res, const acb_t s, const acb_t z, long prec)

    Computes the complementary error function
    `\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)`.
    This function avoids catastrophic cancellation for large positive *z*.

.. function:: void acb_hypgeom_erfi(acb_t res, const acb_t s, const acb_t z, long prec)

    Computes the imaginary error function
    `\operatorname{erfi}(z) = -i\operatorname{erf}(iz)`. This is a trivial wrapper
    of :func:`acb_hypgeom_erf`.

Bessel functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, long prec)

    Computes the Bessel function of the first kind
    via :func:`acb_hypgeom_u_asymp`.
    For all complex `\nu, z`, we have

    .. math ::

        J_{\nu}(z) = \frac{z^{\nu}}{2^{\nu} e^{iz} \Gamma(\nu+1)}
            {}_1F_1(\nu+\tfrac{1}{2}, 2\nu+1, 2iz) = A_{+} B_{+} + A_{-} B_{-}

    where

    .. math ::

        A_{\pm} = z^{\nu} (z^2)^{-\tfrac{1}{2}-\nu} (\mp i z)^{\tfrac{1}{2}+\nu} (2 \pi)^{-1/2} = (\pm iz)^{-1/2-\nu} z^{\nu} (2 \pi)^{-1/2}

        B_{\pm} = e^{\pm i z} U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, \mp 2iz).

    Nicer representations of the factors `A_{\pm}` can be given depending conditionally
    on the parameters. If `\nu + \tfrac{1}{2} = n \in \mathbb{Z}`, we have
    `A_{\pm} = (\pm i)^{n} (2 \pi z)^{-1/2}`.
    And if `\operatorname{Re}(z) > 0`, we have `A_{\pm} = \exp(\mp i [(2\nu+1)/4] \pi) (2 \pi z)^{-1/2}`.

.. function:: void acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, long prec)

    Computes the Bessel function of the first kind from

    .. math ::

        J_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
                     {}_0F_1\left(\nu+1, -\frac{z^2}{4}\right).

.. function:: void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, long prec)

    Computes the Bessel function of the first kind `J_{\nu}(z)` using
    an automatic algorithm choice.

Modified Bessel functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_bessel_k_asymp(acb_t res, const acb_t nu, const acb_t z, long prec)

    Computes the modified Bessel function of the second kind via
    via :func:`acb_hypgeom_u_asymp`. For all `\nu` and all `z \ne 0`, we have

    .. math ::

        K_{\nu}(z) = \left(\frac{\pi}{2z}\right)^{1/2} e^{-z}
            U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, 2z).

.. function:: void acb_hypgeom_bessel_k_0f1_series(acb_poly_t res, const acb_poly_t nu, const acb_poly_t z, long len, long prec)

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

.. function:: void acb_hypgeom_bessel_k_0f1(acb_t res, const acb_t nu, const acb_t z, long prec)

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

.. function:: void acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, long prec)

    Computes the modified Bessel function of the second kind `K_{\nu}(z)` using
    an automatic algorithm choice.

Incomplete gamma functions
-------------------------------------------------------------------------------

.. function:: void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int modified, long prec)

.. function:: void acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s, const acb_t z, int modified, long prec)

.. function:: void acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s, const acb_t z, int modified, long prec)

.. function:: void acb_hypgeom_gamma_upper_singular(acb_t res, long s, const acb_t z, int modified, long prec)

.. function:: void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int modified, long prec)

    Computes the upper incomplete gamma function respectively using

    .. math ::

        \Gamma(s,z) = e^{-z} U(1-s,1-s,z)

        \Gamma(s,z) = \Gamma(s) - \frac{z^s}{s} {}_1F_1(s, s+1, -z)

        \Gamma(s,z) = \Gamma(s) - \frac{z^s e^{-z}}{s} {}_1F_1(1, s+1, z)

        \Gamma(s,z) = \frac{(-1)^n}{n!} (\psi(n+1) - \log(z))
                    + \frac{(-1)^n}{(n+1)!} z \, {}_2F_2(1,1,2,2+n,-z)
                    - z^{-n} \sum_{k=0}^{n-1} \frac{(-z)^k}{(k-n) k!},
                    \quad n = -s \in \mathbb{Z}_{\ge 0}

    and an automatic algorithm choice. The automatic version also handles
    other special input such as `z = 0` and `s = 1, 2, 3`.
    The *singular* version evaluates the finite sum directly and therefore
    assumes that *s* is not too large.
    If *modified* is set, computes the exponential integral
    `z^{-s} \Gamma(s,z) = E_{1-s}(z)` instead.

Exponential and trigonometric integrals
-------------------------------------------------------------------------------

The branch cut conventions of the following functions match Mathematica.

.. function:: void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, long prec)

    Computes the generalized exponential integral `E_s(z)`. This is a
    trivial wrapper of :func:`acb_hypgeom_gamma_upper`.

.. function:: void acb_hypgeom_ei_asymp(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_ei_2f2(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_ei(acb_t res, const acb_t z, long prec)

    Computes the exponential integral `\operatorname{Ei}(z)`, respectively
    using

    .. math ::

        \operatorname{Ei}(z) = -e^z U(1,1,-z) - \log(-z)
            + \frac{1}{2} \left(\log(z) - \log\left(\frac{1}{z}\right) \right)

        \operatorname{Ei}(z) = z {}_2F_2(1, 1; 2, 2; z) + \gamma
            + \frac{1}{2} \left(\log(z) - \log\left(\frac{1}{z}\right) \right)

    and an automatic algorithm choice.

.. function:: void acb_hypgeom_si_asymp(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_si_1f2(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_si(acb_t res, const acb_t z, long prec)

    Computes the sine integral `\operatorname{Si}(z)`, respectively
    using

    .. math ::

        \operatorname{Si}(z) = \frac{i}{2} \left[
            e^{iz} U(1,1,-iz) - e^{-iz} U(1,1,iz) + 
            \log(-iz) - \log(iz) \right]

        \operatorname{Si}(z) = z {}_1F_2(\tfrac{1}{2}; \tfrac{3}{2}, \tfrac{3}{2}; -\tfrac{z^2}{4})

    and an automatic algorithm choice.

.. function:: void acb_hypgeom_ci_asymp(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_ci_2f3(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_ci(acb_t res, const acb_t z, long prec)

    Computes the cosine integral `\operatorname{Ci}(z)`, respectively
    using

    .. math ::

        \operatorname{Ci}(z) = \log(z) - \frac{1}{2} \left[
            e^{iz} U(1,1,-iz) + e^{-iz} U(1,1,iz) + 
            \log(-iz) + \log(iz) \right]

        \operatorname{Ci}(z) = -\tfrac{z^2}{4}
            {}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; -\tfrac{z^2}{4})
            + \log(z) + \gamma

    and an automatic algorithm choice.

.. function:: void acb_hypgeom_shi(acb_t res, const acb_t z, long prec)

    Computes the hyperbolic sine integral
    `\operatorname{Shi}(z) = -i \operatorname{Si}(iz)`.
    This is a trivial wrapper of :func:`acb_hypgeom_si`.

.. function:: void acb_hypgeom_chi_asymp(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_chi_2f3(acb_t res, const acb_t z, long prec)

.. function:: void acb_hypgeom_chi(acb_t res, const acb_t z, long prec)

    Computes the hyperbolic cosine integral `\operatorname{Chi}(z)`, respectively
    using

    .. math ::

        \operatorname{Chi}(z) = -\frac{1}{2} \left[
            e^{z} U(1,1,-z) + e^{-z} U(1,1,z) + 
            \log(-z) - \log(z) \right]

        \operatorname{Chi}(z) = \tfrac{z^2}{4}
            {}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; \tfrac{z^2}{4})
            + \log(z) + \gamma

    and an automatic algorithm choice.

.. function:: void acb_hypgeom_li(acb_t res, const acb_t z, int offset, long prec)

    If *offset* is zero, computes the logarithmic integral
    `\operatorname{li}(z) = \operatorname{Ei}(\log(z))`.

    If *offset* is nonzero, computes the offset logarithmic integral
    `\operatorname{Li}(z) = \operatorname{li}(z) - \operatorname{li}(2)`.

