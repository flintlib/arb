.. _mag:

**mag.h** -- fixed-precision unsigned floating-point numbers for bounds
===============================================================================

The :type:`mag_t` type is an unsigned floating-point type with a
fixed-precision mantissa (30 bits) and unlimited exponent range, suited for
representing and rigorously manipulating magnitude bounds efficiently.
Operations always produce a strict upper or lower bound, but for performance
reasons, no attempt is made to compute the best possible bound
(in general, a result may a few ulps larger/smaller than the optimal value).
The special values zero and positive infinity are supported (but not NaN).
Applications requiring more flexibility (such as correct rounding, or
higher precision) should use the :type:`arf_t` type instead.

The exponent is represented as an :type:`fmpz`. Special "fast" methods are
provided for manipulating :type:`mag_t` objects under the assumption that
no infinities are present and that exponents stay within a certain range.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: mag_struct

    A :type:`mag_struct` holds a mantissa and an exponent.
    Special values are encoded by the mantissa being set to zero.

.. type:: mag_t

    A :type:`mag_t` is defined as an array of length one of type
    :type:`mag_struct`, permitting a :type:`mag_t` to be passed by reference.

Memory management
-------------------------------------------------------------------------------

.. function:: void mag_init(mag_t x)

    Initializes the variable *x* for use. Its value is set to zero.

.. function:: void mag_clear(mag_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: void mag_init_set(mag_t x, const mag_t y)

.. function:: void mag_swap(mag_t x, mag_t y)

.. function:: void mag_set(mag_t x, const mag_t y)

Special values
-------------------------------------------------------------------------------

.. function:: void mag_zero(mag_t x)

    Sets *x* to zero.

.. function:: void mag_inf(mag_t x)

    Sets *x* to positive infinity.

.. function:: int mag_is_special(const mag_t x)

    Returns nonzero iff *x* is zero or positive infinity.

.. function:: int mag_is_zero(const mag_t x)

    Returns nonzero iff *x* is zero.

.. function:: int mag_is_inf(const mag_t x)

    Returns nonzero iff *x* is positive infinity.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void mag_mul(mag_t z, const mag_t x, const mag_t y)

    Sets `z` to an upper bound for `x + y`.

.. function:: void mag_addmul(mag_t z, const mag_t x, const mag_t y)

    Sets `z` to an upper bound for `z + xy`.

.. function:: void mag_add_2exp_fmpz(mag_t z, const mag_t x, const fmpz_t e)

    Sets `z` to an upper bound for `x + 2^e`.

.. function:: void mag_div(mag_t z, const mag_t x, const mag_t y)

    Sets `z` to an upper bound for `x / y`.

Fast versions
-------------------------------------------------------------------------------

.. function:: void mag_fast_init_set(mag_t x, const mag_t y)

.. function:: void mag_fast_zero(mag_t x)

.. function:: int mag_fast_is_zero(const mag_t x)

.. function:: void mag_fast_init_set_arf(mag_t y, const arf_t x)

.. function:: void mag_fast_mul(mag_t z, const mag_t x, const mag_t y)

.. function:: void mag_fast_addmul(mag_t z, const mag_t x, const mag_t y)

.. function:: void mag_fast_add_2exp_si(mag_t z, const mag_t x, long e)

Conversions
-------------------------------------------------------------------------------

These functions are intended for debugging purposes: see the :doc:`arf <arf>`
module for other conversion functions.

.. function:: void mag_set_fmpr(mag_t y, const fmpr_t x)

    Sets *y* to an upper bound for *x*.

.. function:: void mag_get_fmpr(fmpr_t y, const mag_t x)

    Sets *y* to exactly *x*.

