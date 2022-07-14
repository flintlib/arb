.. _acf:

**acf.h** -- complex floating-point numbers
===============================================================================


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: acf_struct

.. type:: acf_t

    An *acf_struct* consists of a pair of *arf_struct*:s.
    An *acf_t* is defined as an array of length one of type
    *acf_struct*, permitting an *acf_t* to be passed by
    reference.

.. type:: acf_ptr

   Alias for ``acf_struct *``, used for vectors of numbers.

.. type:: acf_srcptr

   Alias for ``const acf_struct *``, used for vectors of numbers
   when passed as constant input to functions.

.. macro:: acf_realref(x)

    Macro returning a pointer to the real part of *x* as an *arf_t*.

.. macro:: acf_imagref(x)

    Macro returning a pointer to the imaginary part of *x* as an *arf_t*.


Memory management
-------------------------------------------------------------------------------

.. function:: void acf_init(acf_t x)

    Initializes the variable *x* for use, and sets its value to zero.

.. function:: void acf_clear(acf_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: void acf_swap(acf_t z, acf_t x)

    Swaps *z* and *x* efficiently.

.. function:: slong acf_allocated_bytes(const acf_t x)

    Returns the total number of bytes heap-allocated internally by this object.
    The count excludes the size of the structure itself. Add
    ``sizeof(acf_struct)`` to get the size of the object as a whole.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: arf_ptr acf_real_ptr(acf_t z)
              arf_ptr acf_imag_ptr(acf_t z)

    Returns a pointer to the real or imaginary part of *z*.

.. function:: void acf_set(acf_t z, const acf_t x)

    Sets *z* to the value *x*.

.. function:: int acf_equal(const acf_t x, const acf_t y)

    Returns whether *x* and *y* are equal.
