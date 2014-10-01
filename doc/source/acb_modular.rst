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

The Dedekind eta function
-------------------------------------------------------------------------------

To be done

Jacobi theta functions
-------------------------------------------------------------------------------

To be done

Eisenstein series
-------------------------------------------------------------------------------

To be done

