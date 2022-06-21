.. _fmpzi:

**fmpzi.h** -- Gaussian integers
===============================================================================

This module allows working with elements of the ring `\mathbb{Z}[i]`.
At present, only a minimal interface is provided.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpzi_struct

.. type:: fmpzi_t

    Contains a pairs of integers representing the real and imaginary
    parts. An *fmpzi_t* is defined as an array
    of length one of type *fmpzi_struct*, permitting an *fmpzi_t* to
    be passed by reference.

.. macro:: fmpzi_realref(x)

    Macro giving a pointer to the real part of *x*.

.. macro:: fmpzi_imagref(x)

    Macro giving a pointer to the imaginary part of *x*.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: void fmpzi_init(fmpzi_t x)

.. function:: void fmpzi_clear(fmpzi_t x)

.. function:: int fmpzi_equal(const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_zero(fmpzi_t x)

.. function:: void fmpzi_one(fmpzi_t x)

.. function:: void fmpzi_set(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_set_si_si(fmpzi_t res, slong a, slong b)

.. function:: void fmpzi_swap(fmpzi_t x, fmpzi_t y)

.. function:: void fmpzi_print(const fmpzi_t x)

.. function:: void fmpzi_randtest(fmpzi_t res, flint_rand_t state, mp_bitcnt_t bits)

Arithmetic
-------------------------------------------------------------------------------

.. function:: void fmpzi_conj(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_neg(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_sqr(fmpzi_t res, const fmpzi_t x)

.. function:: void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)

.. function:: void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp)
