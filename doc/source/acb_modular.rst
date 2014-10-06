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


The Dedekind eta function
-------------------------------------------------------------------------------

To be done

.. function:: void acb_modular_addseq_eta(long * exponents, long * aindex, long * bindex, long num)

    Constructs an addition sequence for the first *num* generalized pentagonal
    numbers (excluding zero), i.e. 1, 2, 5, 7, 12, 15, 22, 26, 35, 40 etc.

Jacobi theta functions
-------------------------------------------------------------------------------

To be done

.. function:: void acb_modular_addseq_theta(long * exponents, long * aindex, long * bindex, long num)

    Constructs an addition sequence for the first *num* squares and triangular
    numbers interleaved (excluding zero), i.e. 1, 2, 4, 6, 9, 12, 16, 20, 25, 30 etc.

Eisenstein series
-------------------------------------------------------------------------------

To be done

