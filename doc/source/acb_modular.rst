.. _acb-modular:

**acb_modular.h** -- evaluation of modular forms in the complex numbers
===============================================================================

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

.. function:: void acb_modular_addseq_theta(long * exponents, long * aindex, long * bindex, long num)

    Constructs an addition sequence for the first *num* squares and triangular
    numbers interleaved (excluding zero), i.e. 1, 2, 4, 6, 9, 12, 16, 20, 25, 30 etc.

.. function:: void acb_modular_theta_1234_sum(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t w, int w_is_unit, const acb_t q, long prec)

    Simultaneously evaluates

    .. math ::

        q^{-1/4} \theta_1 = 2 \sum_{n=0}^{\infty} (-1)^n q^{n(n+1)} \sin((2n+1) \pi z)

        q^{-1/4} \theta_2 = 2 \sum_{n=0}^{\infty} q^{n(n+1)} \cos((2n+1) \pi z)

        \theta_3 = 1 + 2 \sum_{n=1}^{\infty} q^{n^2} \cos(2n \pi z)

        \theta_4 = 1 + 2 \sum_{n=1}^{\infty} (-1)^n q^{n^2} \cos(2n \pi z)

    given `w = \exp(\pi i z)` and `q = \exp(\pi i \tau)`.
    If *w_is_unit* is nonzero, *w* is assumed to lie on the unit circle,
    i.e. *z* is assumed to be real.

    Note that the factor `q^{1/4}` is removed from `\theta_1` and `\theta_2`.
    To get the true theta function values, the user has to multiply
    this factor back. This convention avoids an unnecessary root extraction,
    since the user can compute `q^{1/4} = \exp(\pi i \tau / 4)` followed by
    `q = (q^{1/4})^4`, and in many cases when computing products or quotients
    of theta functions, the factor `q^{1/4}` can be eliminated entirely.

    This function is intended for `|q| \ll 1`. It can be called with any
    `q`, but will return useless intervals if convergence is not rapid.
    For general evaluation of theta functions, the user should only call
    this function after applying a suitable modular transformation.

    We consider the sums together, alternatingly updating `(\theta_1, \theta_2)`
    or `(\theta_3, \theta_4)`. For `k = 0, 1, 2, \ldots`, the powers of `q`
    are `\lfloor (k+2)^2 / 4 \rfloor = 1, 2, 4, 6, 9` etc. and the powers of `w` are
    `\pm (k+2) = \pm 2, \pm 3, \pm 4, \ldots` etc.
    For some integer `N \ge 1`, the summation is stopped just before term
    `k = N`. The error can then be bounded as

    .. math ::

        \frac{2 |q|^E \max(|w|,|w^{-1}|)^{N+2}}{1 - |q|^{\lfloor (N+1)/2 \rfloor + 1} \max(|w|,|w^{-1}|)}

    where `E = \lfloor (N+2)^2 / 4 \rfloor`, assuming that the denominator
    is positive.
    This is simply the bound for a geometric series, with the leading
    factor 2 coming from the fact that we sum both negative and positive
    powers of `w`.

    To actually evaluate the series, when `w \ne 1`, we write the even
    cosine terms as `w^{2n} + w^{-2n}`, the odd cosine terms as
    `w (w^{2n} + w^{-2n-2})`, and the sine terms as `w (w^{2n} - w^{-2n-2})`.
    This way we only need even powers of `w` and `w^{-1}`.
    The implementation is not yet optimized for real `z`, in which case
    further work can be saved.

    This function does not permit aliasing between input and output
    arguments.

The Dedekind eta function
-------------------------------------------------------------------------------

.. function:: void acb_modular_addseq_eta(long * exponents, long * aindex, long * bindex, long num)

    Constructs an addition sequence for the first *num* generalized pentagonal
    numbers (excluding zero), i.e. 1, 2, 5, 7, 12, 15, 22, 26, 35, 40 etc.

.. function:: void acb_modular_eta_sum(acb_t eta, const acb_t q, long prec)

    Evaluates the series expansion of the Dedekind eta function
    without the leading 24th root, i.e.

    .. math :: \exp(-\pi i \tau/24) \eta(\tau) = \sum_{n=-\infty}^{\infty} (-1)^n q^{(3n^2-n)/2}

    given `q = \exp(2 \pi i \tau)`.

    This function is intended for `|q| \ll 1`. It can be called with any
    `q`, but will return useless intervals if convergence is not rapid.
    For general evaluation of the eta functions, the user should only call
    this function after applying a suitable modular transformation.

.. function:: acb_modular_epsilon_arg(fmpq_t t, const psl2z_t g)

    Given `g = (a, b; c, d)`, computes a rational number `t` such that
    `\varepsilon(a,b,c,d) = \exp(\pi i t)` is the root of unity in
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
    `j(\tau) = 32 (\theta_2^8 + \theta_3^8 + \theta_4^8)^3 / (\theta_2 \theta_3 \theta_4)^8`
    where `\theta_k` is the respective theta constant evaluated at `\tau`.


