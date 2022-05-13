.. _arb-defs:

**arb-defs.h** -- global definitions for the Arb library
================================================================================

The header ``arb-defs.h`` provides macros and functions used throughout the
Arb-library that are not specific for any Arb-types.


Types, macros and constants
--------------------------------------------------------------------------------

.. macro:: __ARB_VERSION

.. macro:: __ARB_VERSION_MINOR

.. macro:: __ARB_VERSION_PATCHLEVEL

   A macro for the major, minor and patch version of Arb, respectively,
   represented as an integer. In other words, it gives you the first, second or
   third number in the version number.

.. macro:: ARB_VERSION

   A macro for the version of Arb, represented as a string.

.. var:: const char * arb_version

   A constant variable equivalence to :macro:`ARB_VERSION`.

.. macro:: __ARB_RELEASE

   A macro aimed for easily parsing the version number of Arb. It is represented
   as an integer defined by

   .. math ::
      \mathtt{\_\_ARB\_RELEASE}
   =
      10000 \cdot \mathtt{\_\_ARB\_VERSION}
      +
      100 \cdot \mathtt{\_\_ARB\_VERSION\_MINOR}
      +
      \mathtt{\_\_ARB\_VERSION\_PATCHLEVEL}

.. macro:: LIMB_ONE

   The limb with binary representation `(0 0 \cdots 0 1)_{2}`.

.. macro:: LIMB_ONES

   The limb with binary representation `(1 1 \cdots 1 1)_{2}`.

.. macro:: LIMB_TOP

   The limb with binary representation `(1 0 \cdots 0 0)_{2}`.

.. macro:: MASK_LIMB(n, c)

   A macro for removing the `c` lower ones in the bit representation of `n`.

.. macro:: UI_ABS_SI(x)

   Returns the absolute value of `x`.

.. macro:: nn_mul_2x1(r2, r1, r0, a1, a0, b0)

   Given `a = \{a_0, a_1\}` and `b = \{b_0\}`, it sets `r = \{r_0, r_1, r_2\}`
   to the product `a \cdot b`.

.. macro:: nn_mul_2x2(r3, r2, r1, r0, a1, a0, b1, b0)

   Given `a = \{a_0, a_1\}` and `b = \{b_0, b_1\}`, it sets
   `r = \{r_0, r_1, r_2, r_3\}` to the product `a \cdot b`.


Functions
--------------------------------------------------------------------------------

.. function:: double arb_test_multiplier(void)

   Returns the test multiplier. This function is used to determine how many
   tests are going to be performed during ``make check``.

.. function:: int n_zerobits(mp_limb_t e)

   Returns the number of zero bits in the binary representation of `e` *up to
   its most significant bit*.


..
    vim:spell spelllang=en_us:ts=3:sw=3:tw=80:expandtab:
