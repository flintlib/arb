.. _acb-modular:

**acb_modular.h** -- modular forms in the complex numbers
===============================================================================

This module provides methods for numerical evaluation of modular
forms, Jacobi theta functions, and elliptic functions.

In the context of this module, *tau* or `\tau` always denotes an
element of the complex upper half-plane
`\mathbb{H} = \{z \in \mathbb{C} : \operatorname{Im}(z) > 0\}`.
We also often use the variable `q`, variously defined as `q = e^{2 \pi i \tau}`
(usually in relation to modular forms) or `q = e^{\pi i \tau}` (usually
in relation to theta functions) and satisfying `|q| < 1`.
We will clarify the local meaning of `q` every time such a quantity appears as
a function of `\tau`.

As usual, the numerical functions in this module compute strict error
bounds: if *tau* is represented by an :type:`acb_t` whose content
overlaps with the real line (or lies in the lower half-plane),
and *tau* is passed to a function defined only on `\mathbb{H}`, then
the output will have an infinite radius. The analogous behavior holds for
functions requiring `|q| < 1`.

The modular group
-------------------------------------------------------------------------------

.. type:: psl2z_struct

.. type:: psl2z_t

    Represents an element of the modular group `\text{PSL}(2, \mathbb{Z})`,
    namely an integer matrix

    .. math ::

        \begin{pmatrix} a & b \\ c & d \end{pmatrix}

    with `ad-bc = 1`, and with signs canonicalized such that `c \ge 0`,
    and `d > 0` if `c = 0`.
    The struct members *a*, *b*, *c*, *d* are of type *fmpz*.

.. function:: void psl2z_init(psl2z_t g)

    Initializes *g* and set it to the identity element.

.. function:: void psl2z_clear(psl2z_t g)

    Clears *g*.

.. function:: void psl2z_swap(psl2z_t f, psl2z_t g)

    Swaps *f* and *g* efficiently.

.. function:: void psl2z_set(psl2z_t f, const psl2z_t g)

    Sets *f* to a copy of *g*.

.. function:: void psl2z_one(psl2z_t g)

    Sets *g* to the identity element.

.. function:: int psl2z_is_one(const psl2z_t g)

    Returns nonzero iff *g* is the identity element.

.. function:: void psl2z_print(const psl2z_t g)

    Prints *g* to standard output.

.. function:: int psl2z_equal(const psl2z_t f, const psl2z_t g)

    Returns nonzero iff *f* and *g* are equal.

.. function:: void psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g)

    Sets *h* to the product of *f* and *g*, namely the matrix product
    with the signs canonicalized.

.. function:: void psl2z_inv(psl2z_t h, const psl2z_t g)

    Sets *h* to the inverse of *g*.

.. function:: int psl2z_is_correct(const psl2z_t g)

    Returns nonzero iff *g* contains correct data, i.e.
    satisfying `ad-bc = 1`, `c \ge 0`, and `d > 0` if `c = 0`.

.. function:: void psl2z_randtest(psl2z_t g, flint_rand_t state, long bits)

    Sets *g* to a random element of `\text{PSL}(2, \mathbb{Z})`
    with entries of bit length at most *bits*
    (or 1, if *bits* is not positive). We first generate *a* and *d*, compute
    their Bezout coefficients, divide by the GCD, and then correct the signs.

Modular transformations
-------------------------------------------------------------------------------

.. function:: void acb_modular_transform(acb_t w, const psl2z_t g, const acb_t z, long prec)

    Applies the modular transformation *g* to the complex number *z*,
    evaluating

    .. math ::

        w = g z = \frac{az+b}{cz+d}.

.. function:: void acb_modular_fundamental_domain_approx_d(psl2z_t g, double x, double y, double one_minus_eps)

.. function:: void acb_modular_fundamental_domain_approx_arf(psl2z_t g, const arf_t x, const arf_t y, const arf_t one_minus_eps, long prec)

    Attempts to determine a modular transformation *g* that maps the
    complex number `x+yi` to the fundamental domain or just
    slightly outside the fundamental domain, where the target tolerance
    (not a strict bound) is specified by *one_minus_eps*.

    The inputs are assumed to be finite numbers, with *y* positive.

    Uses floating-point iteration, repeatedly applying either
    the transformation `z \gets z + b` or `z \gets -1/z`. The iteration is
    terminated if `|x| \le 1/2` and `x^2 + y^2 \ge 1 - \varepsilon` where
    `1 - \varepsilon` is passed as *one_minus_eps*. It is also terminated
    if too many steps have been taken without convergence, or if the numbers
    end up too large or too small for the working precision.

    The algorithm can fail to produce a satisfactory transformation.
    The output *g* is always set to *some* correct modular transformation,
    but it is up to the user to verify a posteriori that *g* maps `x+yi`
    close enough to the fundamental domain.

.. function:: void acb_modular_fundamental_domain_approx(acb_t w, psl2z_t g, const acb_t z, const arf_t one_minus_eps, long prec)

    Attempts to determine a modular transformation *g* that maps the
    complex number `z` to the fundamental domain or just
    slightly outside the fundamental domain, where the target tolerance
    (not a strict bound) is specified by *one_minus_eps*. It also computes
    the transformed value `w = gz`.

    This function first tries to use
    :func:`acb_modular_fundamental_domain_approx_d` and checks if the
    result is acceptable. If this fails, it calls
    :func:`acb_modular_fundamental_domain_approx_arf` with higher precision.
    Finally, `w = gz` is evaluated by a single application of *g*.

    The algorithm can fail to produce a satisfactory transformation.
    The output *g* is always set to *some* correct modular transformation,
    but it is up to the user to verify a posteriori that `w` is close enough
    to the fundamental domain.

.. function:: int acb_modular_is_in_fundamental_domain(const acb_t z, const arf_t tol, long prec)

    Returns nonzero if it is certainly true that `|z| \ge 1 - \varepsilon` and 
    `|\operatorname{Re}(z)| \le 1/2 + \varepsilon` where `\varepsilon` is
    specified by *tol*. Returns zero if this is false or cannot be determined.

Jacobi theta functions
-------------------------------------------------------------------------------

Unfortunately, there are many inconsistent notational variations for
Jacobi theta functions in the literature. Unless otherwise noted,
we use the functions

.. math ::

    \theta_1(z,\tau) = -i \sum_{n=-\infty}^{\infty} (-1)^n \exp(\pi i [(n + 1/2)^2 \tau + (2n + 1) z])
                     = 2 q_{1/4} \sum_{n=0}^{\infty} (-1)^n q^{n(n+1)} \sin((2n+1) \pi z)

    \theta_2(z,\tau) = \sum_{n=-\infty}^{\infty} \exp(\pi i [(n + 1/2)^2 \tau + (2n + 1) z])
                     = 2 q_{1/4} \sum_{n=0}^{\infty} q^{n(n+1)} \cos((2n+1) \pi z)

    \theta_3(z,\tau) = \sum_{n=-\infty}^{\infty} \exp(\pi i [n^2 \tau + 2n z])
                     = 1 + 2 \sum_{n=1}^{\infty} q^{n^2} \cos(2n \pi z)

    \theta_4(z,\tau) = \sum_{n=-\infty}^{\infty} (-1)^n \exp(\pi i [n^2 \tau + 2n z])
                     = 1 + 2 \sum_{n=1}^{\infty} (-1)^n q^{n^2} \cos(2n \pi z)

where `q = \exp(\pi i \tau)` and `q_{1/4} = \exp(\pi i \tau / 4)`.
Note that many authors write `q_{1/4}` as `q^{1/4}`,
but the principal fourth root `(q)^{1/4} = \exp(\frac{1}{4} \log q)`
differs from `q_{1/4}` in general and some formulas are
only correct if one reads "`q^{1/4} = \exp(\pi i \tau / 4)`".
To avoid confusion, we only write `q^k` when `k` is an integer.

.. function:: void acb_modular_theta_transform(int * R, int * S, int * C, const psl2z_t g)

    We wish to write a theta function with quasiperiod `\tau` in terms
    of a theta function with quasiperiod `\tau' = g \tau`, given
    some `g = (a, b; c, d) \in \text{PSL}(2, \mathbb{Z})`.
    For `i = 0, 1, 2, 3`, this function computes integers `R_i` and `S_i`
    (*R* and *S* should be arrays of length 4)
    and `C \in \{0, 1\}` such that

    .. math ::

        \theta_{1+i}(z,\tau) = \exp(\pi i R_i / 4) \cdot A \cdot B \cdot \theta_{1+S_i}(z',\tau')

    where `z' = z, A = B = 1` if `C = 0`, and

    .. math ::

        z' = \frac{-z}{c \tau + d}, \quad
        A = \sqrt{\frac{i}{c \tau + d}}, \quad
        B = \exp\left(-\pi i c \frac{z^2}{c \tau + d}\right)

    if `C = 1`. Note that `A` is well-defined with the principal branch
    of the square root since `A^2 = i/(c \tau + d)` lies in the right half-plane.

    Firstly, if `c = 0`, we have
    `\theta_i(z, \tau) = \exp(-\pi i b / 4) \theta_i(z, \tau+b)`
    for `i = 1, 2`, whereas
    `\theta_3` and `\theta_4` remain unchanged when `b` is even
    and swap places with each other when `b` is odd.
    In this case we set `C = 0`.

    For an arbitrary `g` with `c > 0`, we set `C = 1`. The general
    transformations are given by Rademacher [Rad1973]_.
    We need the function `\theta_{m,n}(z,\tau)` defined for `m, n \in \mathbb{Z}` by
    (beware of the typos in [Rad1973]_)

    .. math ::

        \theta_{0,0}(z,\tau) = \theta_3(z,\tau), \quad
        \theta_{0,1}(z,\tau) = \theta_4(z,\tau)

        \theta_{1,0}(z,\tau) = \theta_2(z,\tau), \quad
        \theta_{1,1}(z,\tau) = i \theta_1(z,\tau)

        \theta_{m+2,n}(z,\tau) = (-1)^n \theta_{m,n}(z,\tau)

        \theta_{m,n+2}(z,\tau) = \theta_{m,n}(z,\tau).

    Then we may write

    .. math ::

        \theta_1(z,\tau) = \varepsilon_1 A B \theta_1(z', \tau')

        \theta_2(z,\tau) = \varepsilon_2 A B \theta_{1-c,1+a}(z', \tau')

        \theta_3(z,\tau) = \varepsilon_3 A B \theta_{1+d-c,1-b+a}(z', \tau')

        \theta_4(z,\tau) = \varepsilon_4 A B \theta_{1+d,1-b}(z', \tau')

    where `\varepsilon_i` is an 8th root of unity.
    Specifically, if we denote the 24th root of unity
    in the transformation formula of the Dedekind eta
    function by `\varepsilon(a,b,c,d) = \exp(\pi i R(a,b,c,d) / 12)`
    (see :func:`acb_modular_epsilon_arg`), then:

    .. math ::

        \varepsilon_1(a,b,c,d) = \exp(\pi i [R(-d,b,c,-a) + 1] / 4)

        \varepsilon_2(a,b,c,d) = \exp(\pi i [-R(a,b,c,d) + (5+(2-c)a)] / 4)

        \varepsilon_3(a,b,c,d) = \exp(\pi i [-R(a,b,c,d) + (4+(c-d-2)(b-a))] / 4)

        \varepsilon_4(a,b,c,d) = \exp(\pi i [-R(a,b,c,d) + (3-(2+d)b)] / 4)

    These formulas are easily derived from the formulas in [Rad1973]_
    (Rademacher has the transformed/untransformed variables exchanged,
    and his "`\varepsilon`" differs from ours by a constant
    offset in the phase).

.. function:: void acb_modular_addseq_theta(long * exponents, long * aindex, long * bindex, long num)

    Constructs an addition sequence for the first *num* squares and triangular
    numbers interleaved (excluding zero), i.e. 1, 2, 4, 6, 9, 12, 16, 20, 25, 30 etc.

.. function:: void acb_modular_theta_sum(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t w, int w_is_unit, const acb_t q, long len, long prec)

    Simultaneously computes the first *len* coefficients of each of the
    formal power series

    .. math ::

        \theta_1(z+x,\tau) / q_{1/4} \in \mathbb{C}[[x]]

        \theta_2(z+x,\tau) / q_{1/4} \in \mathbb{C}[[x]]

        \theta_3(z+x,\tau) \in \mathbb{C}[[x]]

        \theta_4(z+x,\tau) \in \mathbb{C}[[x]]

    given `w = \exp(\pi i z)` and `q = \exp(\pi i \tau)`, by summing
    a finite truncation of the respective theta function series.
    In particular, with *len* equal to 1, computes the respective
    value of the theta function at the point *z*.
    We require *len* to be positive.
    If *w_is_unit* is nonzero, *w* is assumed to lie on the unit circle,
    i.e. *z* is assumed to be real.

    Note that the factor `q_{1/4}` is removed from `\theta_1` and `\theta_2`.
    To get the true theta function values, the user has to multiply
    this factor back. This convention avoids unnecessary computations,
    since the user can compute `q_{1/4} = \exp(\pi i \tau / 4)` followed by
    `q = (q_{1/4})^4`, and in many cases when computing products or quotients
    of theta functions, the factor `q_{1/4}` can be eliminated entirely.

    This function is intended for `|q| \ll 1`. It can be called with any
    `q`, but will return useless intervals if convergence is not rapid.
    For general evaluation of theta functions, the user should only call
    this function after applying a suitable modular transformation.

    We consider the sums together, alternatingly updating `(\theta_1, \theta_2)`
    or `(\theta_3, \theta_4)`. For `k = 0, 1, 2, \ldots`, the powers of `q`
    are `\lfloor (k+2)^2 / 4 \rfloor = 1, 2, 4, 6, 9` etc. and the powers of `w` are
    `\pm (k+2) = \pm 2, \pm 3, \pm 4, \ldots` etc. The scheme
    is illustrated by the following table:

    .. math ::

        \begin{array}{llll}
               & \theta_1, \theta_2 & q^0 & (w^1 \pm w^{-1}) \\
        k = 0  & \theta_3, \theta_4 & q^1 & (w^2 \pm w^{-2}) \\
        k = 1  & \theta_1, \theta_2 & q^2 & (w^3 \pm w^{-3}) \\
        k = 2  & \theta_3, \theta_4 & q^4 & (w^4 \pm w^{-4}) \\
        k = 3  & \theta_1, \theta_2 & q^6 & (w^5 \pm w^{-5}) \\
        k = 4  & \theta_3, \theta_4 & q^9 & (w^6 \pm w^{-6}) \\
        k = 5  & \theta_1, \theta_2 & q^{12} & (w^7 \pm w^{-7}) \\
        \end{array}

    For some integer `N \ge 1`, the summation is stopped just before term
    `k = N`. Let `Q = |q|`, `W = \max(|w|,|w^{-1}|)`,
    `E = \lfloor (N+2)^2 / 4 \rfloor` and 
    `F = \lfloor (N+1)/2 \rfloor + 1`. The error of the
    zeroth derivative can be bounded as

    .. math ::

        2 Q^E W^{N+2} \left[ 1 + Q^F W + Q^{2F} W^2 + \ldots \right]
        = \frac{2 Q^E W^{N+2}}{1 - Q^F W}

    provided that the denominator is positive (otherwise we set
    the error bound to infinity).
    When *len* is greater than 1, consider the derivative of order *r*.
    The term of index *k* and order *r* picks up a factor of magnitude
    `(k+2)^r` from differentiation of `w^{k+2}` (it also picks up a factor
    `\pi^r`, but we omit this until we rescale the coefficients
    at the end of the computation). Thus we have the error bound

    .. math ::

        2 Q^E W^{N+2} (N+2)^r \left[ 1 + Q^F W \frac{(N+3)^r}{(N+2)^r} + Q^{2F} W^2 \frac{(N+4)^r}{(N+2)^r} + \ldots \right]

    which by the inequality `(1 + m/(N+2))^r \le \exp(mr/(N+2))`
    can be bounded as

    .. math ::

        \frac{2 Q^E W^{N+2} (N+2)^r}{1 - Q^F W \exp(r/(N+2))},

    again valid when the denominator is positive.

    To actually evaluate the series, we write the even
    cosine terms as `w^{2n} + w^{-2n}`, the odd cosine terms as
    `w (w^{2n} + w^{-2n-2})`, and the sine terms as `w (w^{2n} - w^{-2n-2})`.
    This way we only need even powers of `w` and `w^{-1}`.
    The implementation is not yet optimized for real `z`, in which case
    further work can be saved.

    This function does not permit aliasing between input and output
    arguments.

.. function:: void acb_modular_theta_notransform(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, long prec)

    Evaluates the Jacobi theta functions `\theta_i(z,\tau)`, `i = 1, 2, 3, 4`
    simultaneously. This function does not move `\tau` to the fundamental domain.
    This is generally worse than :func:`acb_modular_theta`, but can
    be slightly better for moderate input.

.. function:: void acb_modular_theta(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, long prec)

    Evaluates the Jacobi theta functions `\theta_i(z,\tau)`, `i = 1, 2, 3, 4`
    simultaneously. This function moves `\tau` to the fundamental domain
    before calling :func:`acb_modular_theta_sum`.


The Dedekind eta function
-------------------------------------------------------------------------------

.. function:: void acb_modular_addseq_eta(long * exponents, long * aindex, long * bindex, long num)

    Constructs an addition sequence for the first *num* generalized pentagonal
    numbers (excluding zero), i.e. 1, 2, 5, 7, 12, 15, 22, 26, 35, 40 etc.

.. function:: void acb_modular_eta_sum(acb_t eta, const acb_t q, long prec)

    Evaluates the Dedekind eta function
    without the leading 24th root, i.e.

    .. math :: \exp(-\pi i \tau/12) \eta(\tau) = \sum_{n=-\infty}^{\infty} (-1)^n q^{(3n^2-n)/2}

    given `q = \exp(2 \pi i \tau)`, by summing the defining series.

    This function is intended for `|q| \ll 1`. It can be called with any
    `q`, but will return useless intervals if convergence is not rapid.
    For general evaluation of the eta function, the user should only call
    this function after applying a suitable modular transformation.

.. function:: int acb_modular_epsilon_arg(const psl2z_t g)

    Given `g = (a, b; c, d)`, computes an integer `R` such that
    `\varepsilon(a,b,c,d) = \exp(\pi i R / 12)` is the 24th root of unity in
    the transformation formula for the Dedekind eta function,

    .. math ::

        \eta\left(\frac{a\tau+b}{c\tau+d}\right) = \varepsilon (a,b,c,d)
            \sqrt{c\tau+d} \eta(\tau).

.. function:: void acb_modular_eta(acb_t r, const acb_t tau, long prec)

    Computes the Dedekind eta function `\eta(\tau)` given `\tau` in the upper
    half-plane. This function applies the functional equation to move
    `\tau` to the fundamental domain before calling
    :func:`acb_modular_eta_sum`.

Modular forms
-------------------------------------------------------------------------------

.. function:: void acb_modular_j(acb_t r, const acb_t tau, long prec)

    Computes Klein's j-invariant `j(\tau)` given `\tau` in the upper
    half-plane. The function is normalized so that `j(i) = 1728`.
    We first move `\tau` to the fundamental domain, which does not change
    the value of the function. Then we use the formula
    `j(\tau) = 32 (\theta_2^8+\theta_3^8+\theta_4^8)^3 / (\theta_2 \theta_3 \theta_4)^8` where
    `\theta_i = \theta_i(0,\tau)`.

.. function:: void acb_modular_lambda(acb_t r, const acb_t tau, long prec)

    Computes the lambda function
    `\lambda(\tau) = \theta_2^4(0,\tau) / \theta_3^4(0,\tau)`, which
    is invariant under modular transformations `(a, b; c, d)`
    where `a, d` are odd and `b, c` are even.

.. function:: void acb_modular_delta(acb_t r, const acb_t tau, long prec)

    Computes the modular discriminant `\Delta(\tau) = \eta(\tau)^{24}`,
    which transforms as

    .. math ::

        \Delta\left(\frac{a\tau+b}{c\tau+d}\right) = (c\tau+d)^{12} \Delta(\tau).

    The modular discriminant is sometimes defined with an extra factor
    `(2\pi)^{12}`, which we omit in this implementation.

.. function:: void acb_modular_eisenstein(acb_ptr r, const acb_t tau, long len, long prec)

    Computes simultaneously the first *len* entries in the sequence
    of Eisenstein series `G_4(\tau), G_6(\tau), G_8(\tau), \ldots`,
    defined by

    .. math ::

        G_{2k}(\tau) = \sum_{m^2 + n^2 \ne 0} \frac{1}{(m+n\tau )^{2k}}

    and satisfying

    .. math ::

        G_{2k} \left(\frac{a\tau+b}{c\tau+d}\right) = (c\tau+d)^{2k} G_{2k}(\tau).

    We first evaluate `G_4(\tau)` and `G_6(\tau)` on the fundamental
    domain using theta functions, and then compute the Eisenstein series
    of higher index using a recurrence relation.


Elliptic functions
-------------------------------------------------------------------------------

.. function:: void acb_modular_elliptic_p(acb_t wp, const acb_t z, const acb_t tau, long prec)

    Computes Weierstrass's elliptic function

    .. math ::

        \wp(z, \tau) = \frac{1}{z^2} + \sum_{n^2+m^2 \ne 0}
            \left[ \frac{1}{(z+m+n\tau)^2} - \frac{1}{(m+n\tau)^2} \right]

    which satisfies `\wp(z, \tau) = \wp(z + 1, \tau) = \wp(z + \tau, \tau)`.
    To evaluate the function efficiently, we use the formula

    .. math ::

        \wp(z, \tau) = \pi^2 \theta_2^2(0,\tau) \theta_3^2(0,\tau)
            \frac{\theta_4^2(z,\tau)}{\theta_1^2(z,\tau)} -
            \frac{\pi^2}{3} \left[ \theta_3^4(0,\tau) + \theta_3^4(0,\tau)\right].

.. function:: void acb_modular_elliptic_p_zpx(acb_ptr wp, const acb_t z, const acb_t tau, long len, long prec)

    Computes the formal power series `\wp(z + x, \tau) \in \mathbb{C}[[x]]`,
    truncated to length *len*. In particular, with *len* = 2, simultaneously
    computes `\wp(z, \tau), \wp'(z, \tau)` which together generate
    the field of elliptic functions with periods 1 and `\tau`.

Elliptic integrals
-------------------------------------------------------------------------------

.. function:: void acb_modular_elliptic_k(acb_t w, const acb_t m, long prec)

    Computes the complete elliptic integral of the first kind `K(m)`,
    using the arithmetic-geometric mean: `K(m) = \pi / (2 M(\sqrt{1-m}))`.

.. function:: void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, long len, long prec)

    Sets the coefficients in the array *w* to the power series expansion of the
    complete elliptic integral of the first kind at the point *m* truncated to
    length *len*, i.e. `K(m+x) \in \mathbb{C}[[x]]`.

.. function:: void acb_modular_elliptic_e(acb_t w, const acb_t m, long prec)

    Computes the complete elliptic integral of the second kind `E(m)`,
    which is given by `E(m) = (1-m)(2m K'(m) + K(m))` (where the prime
    denotes a derivative, not a complementary integral).

