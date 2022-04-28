.. _arb:

**arb.h** -- real numbers
===============================================================================

An :type:`arb_t` represents a ball over the real numbers,
that is, an interval `[m \pm r] \equiv [m-r, m+r]` where the midpoint `m` and the
radius `r` are (extended) real numbers and `r` is nonnegative (possibly infinite).
The result of an (approximate) operation done on :type:`arb_t` variables
is a ball which contains the result of the (mathematically exact) operation
applied to any choice of points in the input balls.
In general, the output ball is not the smallest possible.

The precision parameter passed to each function roughly indicates the
precision to which calculations on the midpoint are carried out
(operations on the radius are always done using a fixed, small
precision.)

For arithmetic operations, the precision parameter currently simply
specifies the precision of the corresponding :type:`arf_t` operation.
In the future, the arithmetic might be made faster by incorporating
sloppy rounding (typically equivalent to a loss of 1-2 bits of effective
working precision) when the result is known to be inexact (while still
propagating errors rigorously, of course).
Arithmetic operations done on exact input with exactly
representable output are always guaranteed to produce exact output.

For more complex operations, the precision parameter indicates a minimum
working precision (algorithms might allocate extra internal precision to
attempt to produce an output accurate to the requested number of bits,
especially when the required precision can be estimated easily, but this
is not generally required).

If the precision is increased and the inputs either are exact or are
computed with increased accuracy as well, the output should
converge proportionally, absent any bugs.
The general intended strategy for using ball arithmetic is to add a few
guard bits, and then repeat the calculation as necessary with an
exponentially increasing number of guard bits (Ziv's strategy) until the
result is exact
enough for one's purposes (typically the first attempt will be successful).

The following balls with an infinite or NaN component are permitted,
and may be returned as output from functions.

* The ball `[+\infty \pm c]`, where `c` is finite, represents the point at positive infinity. Such a ball can always be replaced by `[+\infty \pm 0]` while preserving mathematical correctness (this is currently not done automatically by the library).
* The ball `[-\infty \pm c]`, where `c` is finite, represents the point at negative infinity. Such a ball can always be replaced by `[-\infty \pm 0]` while preserving mathematical correctness (this is currently not done automatically by the library).
* The ball `[c \pm \infty]`, where `c` is finite or infinite, represents the whole extended real line `[-\infty,+\infty]`. Such a ball can always be replaced by `[0 \pm \infty]` while preserving mathematical correctness (this is currently not done automatically by the library). Note that there is no way to represent a half-infinite interval such as `[0,\infty]`.
* The ball `[\operatorname{NaN} \pm c]`, where `c` is finite or infinite, represents an indeterminate value (the value could be any extended real number, or it could represent a function being evaluated outside its domain of definition, for example where the result would be complex). Such an indeterminate ball can always be replaced by `[\operatorname{NaN} \pm \infty]` while preserving mathematical correctness (this is currently not done automatically by the library).

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: arb_struct

.. type:: arb_t

    An :type:`arb_struct` consists of an :type:`arf_struct` (the midpoint) and
    a :type:`mag_struct` (the radius).
    An :type:`arb_t` is defined as an array of length one of type
    :type:`arb_struct`, permitting an :type:`arb_t` to be passed by
    reference.

.. type:: arb_ptr

   Alias for ``arb_struct *``, used for vectors of numbers.

.. type:: arb_srcptr

   Alias for ``const arb_struct *``, used for vectors of numbers
   when passed as constant input to functions.

.. macro:: arb_midref(x)

    Macro returning a pointer to the midpoint of *x* as an :type:`arf_t`.

.. macro:: arb_radref(x)

    Macro returning a pointer to the radius of *x* as a :type:`mag_t`.

Memory management
-------------------------------------------------------------------------------

.. function:: void arb_init(arb_t x)

    Initializes the variable *x* for use. Its midpoint and radius are both
    set to zero.

.. function:: void arb_clear(arb_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

.. function:: arb_ptr _arb_vec_init(slong n)

    Returns a pointer to an array of *n* initialized :type:`arb_struct`
    entries.

.. function:: void _arb_vec_clear(arb_ptr v, slong n)

    Clears an array of *n* initialized :type:`arb_struct` entries.

.. function:: void arb_swap(arb_t x, arb_t y)

    Swaps *x* and *y* efficiently.

.. function:: slong arb_allocated_bytes(const arb_t x)

    Returns the total number of bytes heap-allocated internally by this object.
    The count excludes the size of the structure itself. Add
    ``sizeof(arb_struct)`` to get the size of the object as a whole.

.. function:: slong _arb_vec_allocated_bytes(arb_srcptr vec, slong len)

    Returns the total number of bytes allocated for this vector, i.e. the
    space taken up by the vector itself plus the sum of the internal heap
    allocation sizes for all its member elements.

.. function:: double _arb_vec_estimate_allocated_bytes(slong len, slong prec)

    Estimates the number of bytes that need to be allocated for a vector of
    *len* elements with *prec* bits of precision, including the space for
    internal limb data.
    This function returns a *double* to avoid overflow issues when both
    *len* and *prec* are large.

    This is only an approximation of the physical memory that will be used
    by an actual vector. In practice, the space varies with the content
    of the numbers; for example, zeros and small integers require no
    internal heap allocation even if the precision is huge.
    The estimate assumes that exponents will not be bignums.
    The actual amount may also be higher or lower due to overhead in the
    memory allocator or overcommitment by the operating system.

Assignment and rounding
-------------------------------------------------------------------------------

.. function:: void arb_set(arb_t y, const arb_t x)

.. function:: void arb_set_arf(arb_t y, const arf_t x)

.. function:: void arb_set_si(arb_t y, slong x)

.. function:: void arb_set_ui(arb_t y, ulong x)

.. function:: void arb_set_d(arb_t y, double x)

.. function:: void arb_set_fmpz(arb_t y, const fmpz_t x)

    Sets *y* to the value of *x* without rounding.

.. function:: void arb_set_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t e)

    Sets *y* to `x \cdot 2^e`.

.. function:: void arb_set_round(arb_t y, const arb_t x, slong prec)

.. function:: void arb_set_round_fmpz(arb_t y, const fmpz_t x, slong prec)

    Sets *y* to the value of *x*, rounded to *prec* bits in the direction
    towards zero.

.. function:: void arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t e, slong prec)

    Sets *y* to `x \cdot 2^e`, rounded to *prec* bits in the direction
    towards zero.

.. function:: void arb_set_fmpq(arb_t y, const fmpq_t x, slong prec)

    Sets *y* to the rational number *x*, rounded to *prec* bits in the direction
    towards zero.

.. function:: int arb_set_str(arb_t res, const char * inp, slong prec)

    Sets *res* to the value specified by the human-readable string *inp*.
    The input may be a decimal floating-point literal,
    such as "25", "0.001", "7e+141" or "-31.4159e-1", and may also consist
    of two such literals separated by the symbol "+/-" and optionally
    enclosed in brackets, e.g. "[3.25 +/- 0.0001]", or simply
    "[+/- 10]" with an implicit zero midpoint.
    The output is rounded to *prec* bits, and if the binary-to-decimal
    conversion is inexact, the resulting error is added to the radius.

    The symbols "inf" and "nan" are recognized (a nan midpoint results in an
    indeterminate interval, with infinite radius).

    Returns 0 if successful and nonzero if unsuccessful. If unsuccessful,
    the result is set to an indeterminate interval.

.. function:: char * arb_get_str(const arb_t x, slong n, ulong flags)

    Returns a nice human-readable representation of *x*, with at most *n*
    digits of the midpoint printed.

    With default flags, the output can be parsed back with :func:`arb_set_str`,
    and this is guaranteed to produce an interval containing the original
    interval *x*.

    By default, the output is rounded so that the value given for the
    midpoint is correct up to 1 ulp (unit in the last decimal place).

    If *ARB_STR_MORE* is added to *flags*, more (possibly incorrect)
    digits may be printed.

    If *ARB_STR_NO_RADIUS* is added to *flags*, the radius is not
    included in the output. Unless *ARB_STR_MORE* is set, the output is
    rounded so that the midpoint is correct to 1 ulp. As a special case,
    if there are no significant digits after rounding, the result will
    be shown as ``0e+n``, meaning that the result is between
    ``-1e+n`` and ``1e+n`` (following the contract that the output is
    correct to within one unit in the only shown digit).

    By adding a multiple *m* of *ARB_STR_CONDENSE* to *flags*, strings
    of more than three times *m* consecutive digits are condensed, only
    printing the leading and trailing *m* digits along with
    brackets indicating the number of digits omitted
    (useful when computing values to extremely high precision).

Assignment of special values
-------------------------------------------------------------------------------

.. function:: void arb_zero(arb_t x)

    Sets *x* to zero.

.. function:: void arb_one(arb_t f)

    Sets *x* to the exact integer 1.

.. function:: void arb_pos_inf(arb_t x)

    Sets *x* to positive infinity, with a zero radius.

.. function:: void arb_neg_inf(arb_t x)

    Sets *x* to negative infinity, with a zero radius.

.. function:: void arb_zero_pm_inf(arb_t x)

    Sets *x* to `[0 \pm \infty]`, representing the whole extended real line.

.. function:: void arb_indeterminate(arb_t x)

    Sets *x* to `[\operatorname{NaN} \pm \infty]`, representing
    an indeterminate result.

.. function:: void arb_zero_pm_one(arb_t x)

    Sets *x* to the interval `[0 \pm 1]`.

.. function:: void arb_unit_interval(arb_t x)

    Sets *x* to the interval `[0, 1]`.

Input and output
-------------------------------------------------------------------------------

The *arb_print...* functions print to standard output, while
*arb_fprint...* functions print to the stream *file*.

.. function:: void arb_print(const arb_t x)

.. function:: void arb_fprint(FILE * file, const arb_t x)

    Prints the internal representation of *x*.

.. function:: void arb_printd(const arb_t x, slong digits)

.. function:: void arb_fprintd(FILE * file, const arb_t x, slong digits)

    Prints *x* in decimal. The printed value of the radius is not adjusted
    to compensate for the fact that the binary-to-decimal conversion
    of both the midpoint and the radius introduces additional error.

.. function:: void arb_printn(const arb_t x, slong digits, ulong flags)

.. function:: void arb_fprintn(FILE * file, const arb_t x, slong digits, ulong flags)

    Prints a nice decimal representation of *x*.
    By default, the output shows the midpoint with a guaranteed error of at
    most one unit in the last decimal place. In addition, an explicit error
    bound is printed so that the displayed decimal interval is guaranteed to
    enclose *x*.
    See :func:`arb_get_str` for details.

.. function:: char * arb_dump_str(const arb_t x)

    Returns a serialized representation of *x* as a null-terminated
    ASCII string that can be read by :func:`arb_load_str`. The format consists
    of four hexadecimal integers representing the midpoint mantissa,
    midpoint exponent, radius mantissa and radius exponent (with special
    values to indicate zero, infinity and NaN values),
    separated by single spaces. The returned string needs to be deallocated
    with *flint_free*.

.. function:: int arb_load_str(arb_t x, const char * str)

    Sets *x* to the serialized representation given in *str*. Returns a
    nonzero value if *str* is not formatted correctly (see :func:`arb_dump_str`).

.. function:: int arb_dump_file(FILE * stream, const arb_t x)

    Writes a serialized ASCII representation of *x* to *stream* in a form that
    can be read by :func:`arb_load_file`. Returns a nonzero value if the data
    could not be written.

.. function:: int arb_load_file(arb_t x, FILE * stream)

    Reads *x* from a serialized ASCII representation in *stream*. Returns a
    nonzero value if the data is not
    formatted correctly or the read failed. Note that the data is assumed to be
    delimited by a whitespace or end-of-file, i.e., when writing multiple
    values with :func:`arb_dump_file` make sure to insert a whitespace to
    separate consecutive values.

    It is possible to serialize and deserialize a vector as follows
    (warning: without error handling):

    .. code-block:: c

        fp = fopen("data.txt", "w");
        for (i = 0; i < n; i++)
        {
            arb_dump_file(fp, vec + i);
            fprintf(fp, "\n");    // or any whitespace character
        }
        fclose(fp);

        fp = fopen("data.txt", "r");
        for (i = 0; i < n; i++)
        {
            arb_load_file(vec + i, fp);
        }
        fclose(fp);


Random number generation
-------------------------------------------------------------------------------

.. function:: void arb_randtest(arb_t x, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random ball. The midpoint and radius will both be finite.

.. function:: void arb_randtest_exact(arb_t x, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random number with zero radius.

.. function:: void arb_randtest_precise(arb_t x, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random number with radius around `2^{-\text{prec}}`
    the magnitude of the midpoint.

.. function:: void arb_randtest_wide(arb_t x, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random number with midpoint and radius chosen independently,
    possibly giving a very large interval.

.. function:: void arb_randtest_special(arb_t x, flint_rand_t state, slong prec, slong mag_bits)

    Generates a random interval, possibly having NaN or an infinity
    as the midpoint and possibly having an infinite radius.

.. function:: void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, slong bits)

    Sets *q* to a random rational number from the interval represented by *x*.
    A denominator is chosen by multiplying the binary denominator of *x*
    by a random integer up to *bits* bits.

    The outcome is undefined if the midpoint or radius of *x* is non-finite,
    or if the exponent of the midpoint or radius is so large or small
    that representing the endpoints as exact rational numbers would
    cause overflows.

.. function:: void arb_urandom(arb_t x, flint_rand_t state, slong prec, arf_rnd_t rnd)

    Sets *x* to a uniformly distributed random number in the interval
    `[0, 1]`. The method uses rounding from integers to floats, hence the
    radius might not be `0`.

Radius and interval operations
-------------------------------------------------------------------------------

.. function:: void arb_get_mid_arb(arb_t m, const arb_t x)

    Sets *m* to the midpoint of *x*.

.. function:: void arb_get_rad_arb(arb_t r, const arb_t x)

    Sets *r* to the radius of *x*.

.. function:: void arb_add_error_arf(arb_t x, const arf_t err)

.. function:: void arb_add_error_mag(arb_t x, const mag_t err)

.. function:: void arb_add_error(arb_t x, const arb_t err)

    Adds the absolute value of *err* to the radius of *x* (the operation
    is done in-place).

.. function:: void arb_add_error_2exp_si(arb_t x, slong e)

.. function:: void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t e)

    Adds `2^e` to the radius of *x*.

.. function:: void arb_union(arb_t z, const arb_t x, const arb_t y, slong prec)

    Sets *z* to a ball containing both *x* and *y*.

.. function:: int arb_intersection(arb_t z, const arb_t x, const arb_t y, slong prec)

    If *x* and *y* overlap according to :func:`arb_overlaps`,
    then *z* is set to a ball containing the intersection of *x* and *y*
    and a nonzero value is returned.
    Otherwise zero is returned and the value of *z* is undefined.
    If *x* or *y* contains NaN, the result is NaN.

.. function:: void arb_nonnegative_part(arb_t res, const arb_t x)

    Sets *res* to the intersection of *x* with `[0,\infty]`. If *x* is
    nonnegative, an exact copy is made. If *x* is finite and contains negative
    numbers, an interval of the form `[r/2 \pm r/2]` is produced, which
    certainly contains no negative points.
    In the special case when *x* is strictly negative, *res* is set to zero.

.. function:: void arb_get_abs_ubound_arf(arf_t u, const arb_t x, slong prec)

    Sets *u* to the upper bound for the absolute value of *x*,
    rounded up to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void arb_get_abs_lbound_arf(arf_t u, const arb_t x, slong prec)

    Sets *u* to the lower bound for the absolute value of *x*,
    rounded down to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void arb_get_ubound_arf(arf_t u, const arb_t x, slong prec)

    Sets *u* to the upper bound for the value of *x*,
    rounded up to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void arb_get_lbound_arf(arf_t u, const arb_t x, slong prec)

    Sets *u* to the lower bound for the value of *x*,
    rounded down to *prec* bits. If *x* contains NaN, the result is NaN.

.. function:: void arb_get_mag(mag_t z, const arb_t x)

    Sets *z* to an upper bound for the absolute value of *x*. If *x* contains
    NaN, the result is positive infinity.

.. function:: void arb_get_mag_lower(mag_t z, const arb_t x)

    Sets *z* to a lower bound for the absolute value of *x*. If *x* contains
    NaN, the result is zero.

.. function:: void arb_get_mag_lower_nonnegative(mag_t z, const arb_t x)

    Sets *z* to a lower bound for the signed value of *x*, or zero
    if *x* overlaps with the negative half-axis. If *x* contains NaN,
    the result is zero.

.. function:: void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x)

    Computes the exact interval represented by *x*, in the form of an integer
    interval multiplied by a power of two, i.e. `x = [a, b] \times 2^{\text{exp}}`.
    The result is normalized by removing common trailing zeros
    from *a* and *b*.

    This method aborts if *x* is infinite or NaN, or if the difference between
    the exponents of the midpoint and the radius is so large that allocating
    memory for the result fails.

    Warning: this method will allocate a huge amount of memory to store
    the result if the exponent difference is huge. Memory allocation could
    succeed even if the required space is far larger than the physical
    memory available on the machine, resulting in swapping. It is recommended
    to check that the midpoint and radius of *x* both are within a
    reasonable range before calling this method.

.. function:: void arb_set_interval_mag(arb_t x, const mag_t a, const mag_t b, slong prec)

.. function:: void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, slong prec)

.. function:: void arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, slong prec)

    Sets *x* to a ball containing the interval `[a, b]`. We
    require that `a \le b`.

.. function:: void arb_set_interval_neg_pos_mag(arb_t x, const mag_t a, const mag_t b, slong prec)

    Sets *x* to a ball containing the interval `[-a, b]`.

.. function:: void arb_get_interval_arf(arf_t a, arf_t b, const arb_t x, slong prec)

.. function:: void arb_get_interval_mpfr(mpfr_t a, mpfr_t b, const arb_t x)

    Constructs an interval `[a, b]` containing the ball *x*. The MPFR version
    uses the precision of the output variables.

.. function:: slong arb_rel_error_bits(const arb_t x)

    Returns the effective relative error of *x* measured in bits, defined as
    the difference between the position of the top bit in the radius
    and the top bit in the midpoint, plus one.
    The result is clamped between plus/minus *ARF_PREC_EXACT*.

.. function:: slong arb_rel_accuracy_bits(const arb_t x)

    Returns the effective relative accuracy of *x* measured in bits,
    equal to the negative of the return value from :func:`arb_rel_error_bits`.

.. function:: slong arb_rel_one_accuracy_bits(const arb_t x)

    Given a ball with midpoint *m* and radius *r*, returns an approximation of
    the relative accuracy of `[\max(1,|m|) \pm r]` measured in bits.

.. function:: slong arb_bits(const arb_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of the midpoint of *x*, i.e. the minimum precision
    sufficient to represent *x* exactly. Returns 0 if the midpoint
    of *x* is a special value.

.. function:: void arb_trim(arb_t y, const arb_t x)

    Sets *y* to a trimmed copy of *x*: rounds *x* to a number of bits
    equal to the accuracy of *x* (as indicated by its radius),
    plus a few guard bits. The resulting ball is guaranteed to
    contain *x*, but is more economical if *x* has
    less than full accuracy.

.. function:: int arb_get_unique_fmpz(fmpz_t z, const arb_t x)

    If *x* contains a unique integer, sets *z* to that value and returns
    nonzero. Otherwise (if *x* represents no integers or more than one integer),
    returns zero.

    This method aborts if there is a unique integer but that integer
    is so large that allocating memory for the result fails.

    Warning: this method will allocate a huge amount of memory to store
    the result if there is a unique integer and that integer is huge.
    Memory allocation could succeed even if the required space is far
    larger than the physical memory available on the machine, resulting
    in swapping. It is recommended to check that the midpoint of *x* is
    within a reasonable range before calling this method.

.. function:: void arb_floor(arb_t y, const arb_t x, slong prec)

.. function:: void arb_ceil(arb_t y, const arb_t x, slong prec)

    Sets *y* to a ball containing `\lfloor x \rfloor` and `\lceil x \rceil`
    respectively, with the midpoint of *y* rounded to at most *prec* bits.

.. function:: void arb_get_fmpz_mid_rad_10exp(fmpz_t mid, fmpz_t rad, fmpz_t exp, const arb_t x, slong n)

    Assuming that *x* is finite and not exactly zero, computes integers *mid*,
    *rad*, *exp* such that `x \in [m-r, m+r] \times 10^e` and such that the
    larger out of *mid* and *rad* has at least *n* digits plus a few guard
    digits. If *x* is infinite or exactly zero, the outputs are all set
    to zero.

.. function:: int arb_can_round_arf(const arb_t x, slong prec, arf_rnd_t rnd)

.. function:: int arb_can_round_mpfr(const arb_t x, slong prec, mpfr_rnd_t rnd)

    Returns nonzero if rounding the midpoint of *x* to *prec* bits in
    the direction *rnd* is guaranteed to give the unique correctly
    rounded floating-point approximation for the real number represented by *x*.

    In other words, if this function returns nonzero, applying
    :func:`arf_set_round`, or :func:`arf_get_mpfr`, or :func:`arf_get_d`
    to the midpoint of *x* is guaranteed to return a correctly rounded *arf_t*,
    *mpfr_t* (provided that *prec* is the precision of the output variable),
    or *double* (provided that *prec* is 53).
    Moreover, :func:`arf_get_mpfr` is guaranteed to return the correct ternary
    value according to MPFR semantics.

    Note that the *mpfr* version of this function takes an MPFR rounding mode
    symbol as input, while the *arf* version takes an *arf* rounding mode
    symbol. Otherwise, the functions are identical.

    This function may perform a fast, inexact test; that is, it may return
    zero in some cases even when correct rounding actually is possible.

    To be conservative, zero is returned when *x* is non-finite, even if it
    is an "exact" infinity.

Comparisons
-------------------------------------------------------------------------------

.. function:: int arb_is_zero(const arb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both zero.

.. function:: int arb_is_nonzero(const arb_t x)

    Returns nonzero iff zero is not contained in the interval represented
    by *x*.

.. function:: int arb_is_one(const arb_t f)

    Returns nonzero iff *x* is exactly 1.

.. function:: int arb_is_finite(const arb_t x)

    Returns nonzero iff the midpoint and radius of *x* are both finite
    floating-point numbers, i.e. not infinities or NaN.

.. function:: int arb_is_exact(const arb_t x)

    Returns nonzero iff the radius of *x* is zero.

.. function:: int arb_is_int(const arb_t x)

    Returns nonzero iff *x* is an exact integer.

.. function:: int arb_is_int_2exp_si(const arb_t x, slong e)

    Returns nonzero iff *x* exactly equals `n 2^e` for some integer *n*.

.. function:: int arb_equal(const arb_t x, const arb_t y)

    Returns nonzero iff *x* and *y* are equal as balls, i.e. have both the
    same midpoint and radius.

    Note that this is not the same thing as testing whether both
    *x* and *y* certainly represent the same real number, unless
    either *x* or *y* is exact (and neither contains NaN).
    To test whether both operands *might* represent the same mathematical
    quantity, use :func:`arb_overlaps` or :func:`arb_contains`,
    depending on the circumstance.

.. function:: int arb_equal_si(const arb_t x, slong y)

    Returns nonzero iff *x* is equal to the integer *y*.

.. function:: int arb_is_positive(const arb_t x)

.. function:: int arb_is_nonnegative(const arb_t x)

.. function:: int arb_is_negative(const arb_t x)

.. function:: int arb_is_nonpositive(const arb_t x)

    Returns nonzero iff all points *p* in the interval represented by *x*
    satisfy, respectively, `p > 0`, `p \ge 0`, `p < 0`, `p \le 0`.
    If *x* contains NaN, returns zero.

.. function:: int arb_overlaps(const arb_t x, const arb_t y)

    Returns nonzero iff *x* and *y* have some point in common.
    If either *x* or *y* contains NaN, this function always returns nonzero
    (as a NaN could be anything, it could in particular contain any
    number that is included in the other operand).

.. function:: int arb_contains_arf(const arb_t x, const arf_t y)

.. function:: int arb_contains_fmpq(const arb_t x, const fmpq_t y)

.. function:: int arb_contains_fmpz(const arb_t x, const fmpz_t y)

.. function:: int arb_contains_si(const arb_t x, slong y)

.. function:: int arb_contains_mpfr(const arb_t x, const mpfr_t y)

.. function:: int arb_contains(const arb_t x, const arb_t y)

    Returns nonzero iff the given number (or ball) *y* is contained in
    the interval represented by *x*.

    If *x* contains NaN, this function always returns nonzero (as it
    could represent anything, and in particular could represent all
    the points included in *y*).
    If *y* contains NaN and *x* does not, it always returns zero.

.. function:: int arb_contains_int(const arb_t x)

    Returns nonzero iff the interval represented by *x* contains an integer.

.. function:: int arb_contains_zero(const arb_t x)

.. function:: int arb_contains_negative(const arb_t x)

.. function:: int arb_contains_nonpositive(const arb_t x)

.. function:: int arb_contains_positive(const arb_t x)

.. function:: int arb_contains_nonnegative(const arb_t x)

    Returns nonzero iff there is any point *p* in the interval represented
    by *x* satisfying, respectively, `p = 0`, `p < 0`, `p \le 0`, `p > 0`, `p \ge 0`.
    If *x* contains NaN, returns nonzero.

.. function:: int arb_contains_interior(const arb_t x, const arb_t y)

    Tests if *y* is contained in the interior of *x*; that is, contained
    in *x* and not touching either endpoint.

.. function:: int arb_eq(const arb_t x, const arb_t y)

.. function:: int arb_ne(const arb_t x, const arb_t y)

.. function:: int arb_lt(const arb_t x, const arb_t y)

.. function:: int arb_le(const arb_t x, const arb_t y)

.. function:: int arb_gt(const arb_t x, const arb_t y)

.. function:: int arb_ge(const arb_t x, const arb_t y)

    Respectively performs the comparison `x = y`, `x \ne y`,
    `x < y`, `x \le y`, `x > y`, `x \ge y` in a mathematically meaningful way.
    If the comparison `t \, (\operatorname{op}) \, u` holds for all
    `t \in x` and all `u \in y`, returns 1.
    Otherwise, returns 0.

    The balls *x* and *y* are viewed as subintervals of the extended real line.
    Note that balls that are formally different can compare as equal
    under this definition: for example, `[-\infty \pm 3] = [-\infty \pm 0]`.
    Also `[-\infty] \le [\infty \pm \infty]`.

    The output is always 0 if either input has NaN as midpoint.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void arb_neg(arb_t y, const arb_t x)

.. function:: void arb_neg_round(arb_t y, const arb_t x, slong prec)

    Sets *y* to the negation of *x*.

.. function:: void arb_abs(arb_t y, const arb_t x)

    Sets *y* to the absolute value of *x*. No attempt is made to improve the
    interval represented by *x* if it contains zero.

.. function:: void arb_sgn(arb_t y, const arb_t x)

    Sets *y* to the sign function of *x*. The result is `[0 \pm 1]` if
    *x* contains both zero and nonzero numbers.

.. function:: int arb_sgn_nonzero(const arb_t x)

    Returns 1 if *x* is strictly positive, -1 if *x* is strictly negative,
    and 0 if *x* is zero or a ball containing zero so that its sign
    is not determined.

.. function:: void arb_min(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_max(arb_t z, const arb_t x, const arb_t y, slong prec)

    Sets *z* respectively to the minimum and the maximum of *x* and *y*.

.. function:: void arb_add(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_add_arf(arb_t z, const arb_t x, const arf_t y, slong prec)

.. function:: void arb_add_ui(arb_t z, const arb_t x, ulong y, slong prec)

.. function:: void arb_add_si(arb_t z, const arb_t x, slong y, slong prec)

.. function:: void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

    Sets `z = x + y`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t m, const fmpz_t e, slong prec)

    Sets `z = x + m \cdot 2^e`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_sub(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, slong prec)

.. function:: void arb_sub_ui(arb_t z, const arb_t x, ulong y, slong prec)

.. function:: void arb_sub_si(arb_t z, const arb_t x, slong y, slong prec)

.. function:: void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

    Sets `z = x - y`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_mul(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_mul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)

.. function:: void arb_mul_si(arb_t z, const arb_t x, slong y, slong prec)

.. function:: void arb_mul_ui(arb_t z, const arb_t x, ulong y, slong prec)

.. function:: void arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

    Sets `z = x \cdot y`, rounded to *prec* bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_mul_2exp_si(arb_t y, const arb_t x, slong e)

.. function:: void arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)

    Sets *y* to *x* multiplied by `2^e`.

.. function:: void arb_addmul(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)

.. function:: void arb_addmul_si(arb_t z, const arb_t x, slong y, slong prec)

.. function:: void arb_addmul_ui(arb_t z, const arb_t x, ulong y, slong prec)

.. function:: void arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

    Sets `z = z + x \cdot y`, rounded to prec bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_submul(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_submul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)

.. function:: void arb_submul_si(arb_t z, const arb_t x, slong y, slong prec)

.. function:: void arb_submul_ui(arb_t z, const arb_t x, ulong y, slong prec)

.. function:: void arb_submul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

    Sets `z = z - x \cdot y`, rounded to prec bits. The precision can be
    *ARF_PREC_EXACT* provided that the result fits in memory.

.. function:: void arb_fma(arb_t res, const arb_t x, const arb_t y, const arb_t z, slong prec)
              void arb_fma_arf(arb_t res, const arb_t x, const arf_t y, const arb_t z, slong prec)
              void arb_fma_si(arb_t res, const arb_t x, slong y, const arb_t z, slong prec)
              void arb_fma_ui(arb_t res, const arb_t x, ulong y, const arb_t z, slong prec)
              void arb_fma_fmpz(arb_t res, const arb_t x, const fmpz_t y, const arb_t z, slong prec)

    Sets *res* to `x \cdot y + z`. This is equivalent to an *addmul* except
    that *res* and *z* can be separate variables.

.. function:: void arb_inv(arb_t z, const arb_t x, slong prec)

    Sets *z* to `1 / x`.

.. function:: void arb_div(arb_t z, const arb_t x, const arb_t y, slong prec)

.. function:: void arb_div_arf(arb_t z, const arb_t x, const arf_t y, slong prec)

.. function:: void arb_div_si(arb_t z, const arb_t x, slong y, slong prec)

.. function:: void arb_div_ui(arb_t z, const arb_t x, ulong y, slong prec)

.. function:: void arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)

.. function:: void arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, slong prec)

.. function:: void arb_ui_div(arb_t z, ulong x, const arb_t y, slong prec)

    Sets `z = x / y`, rounded to *prec* bits. If *y* contains zero, *z* is
    set to `0 \pm \infty`. Otherwise, error propagation uses the rule

    .. math ::
        \left| \frac{x}{y} - \frac{x+\xi_1 a}{y+\xi_2 b} \right| =
        \left|\frac{x \xi_2 b - y \xi_1 a}{y (y+\xi_2 b)}\right| \le
        \frac{|xb|+|ya|}{|y| (|y|-b)}

    where `-1 \le \xi_1, \xi_2 \le 1`, and
    where the triangle inequality has been applied to the numerator and
    the reverse triangle inequality has been applied to the denominator.

.. function:: void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, slong prec)

    Sets `z = x / (2^n - 1)`, rounded to *prec* bits.

Dot product
-------------------------------------------------------------------------------

.. function:: void arb_dot_precise(arb_t res, const arb_t s, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
              void arb_dot_simple(arb_t res, const arb_t s, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
              void arb_dot(arb_t res, const arb_t s, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)

    Computes the dot product of the vectors *x* and *y*, setting
    *res* to `s + (-1)^{subtract} \sum_{i=0}^{len-1} x_i y_i`.

    The initial term *s* is optional and can be
    omitted by passing *NULL* (equivalently, `s = 0`).
    The parameter *subtract* must be 0 or 1.
    The length *len* is allowed to be negative, which is equivalent
    to a length of zero.
    The parameters *xstep* or *ystep* specify a step length for
    traversing subsequences of the vectors *x* and *y*; either can be
    negative to step in the reverse direction starting from
    the initial pointer.
    Aliasing is allowed between *res* and *s* but not between
    *res* and the entries of *x* and *y*.

    The default version determines the optimal precision for each term
    and performs all internal calculations using mpn arithmetic
    with minimal overhead. This is the preferred way to compute a
    dot product; it is generally much faster and more precise
    than a simple loop.

    The *simple* version performs fused multiply-add operations in
    a simple loop. This can be used for
    testing purposes and is also used as a fallback by the
    default version when the exponents are out of range
    for the optimized code.

    The *precise* version computes the dot product exactly up to the
    final rounding. This can be extremely slow and is only intended
    for testing.

.. function:: void arb_approx_dot(arb_t res, const arb_t s, int subtract, arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)

    Computes an approximate dot product *without error bounds*.
    The radii of the inputs are ignored (only the midpoints are read)
    and only the midpoint of the output is written.

.. function:: void arb_dot_ui(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
              void arb_dot_si(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const slong * y, slong ystep, slong len, slong prec)
              void arb_dot_uiui(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
              void arb_dot_siui(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec)
              void arb_dot_fmpz(arb_t res, const arb_t initial, int subtract, arb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec)

    Equivalent to :func:`arb_dot`, but with integers in the array *y*.
    The *uiui* and *siui* versions take an array of double-limb integers
    as input; the *siui* version assumes that these represent signed
    integers in two's complement form.


Powers and roots
-------------------------------------------------------------------------------

.. function:: void arb_sqrt(arb_t z, const arb_t x, slong prec)

.. function:: void arb_sqrt_arf(arb_t z, const arf_t x, slong prec)

.. function:: void arb_sqrt_fmpz(arb_t z, const fmpz_t x, slong prec)

.. function:: void arb_sqrt_ui(arb_t z, ulong x, slong prec)

    Sets *z* to the square root of *x*, rounded to *prec* bits.

    If `x = m \pm x` where `m \ge r \ge 0`, the propagated error is bounded by
    `\sqrt{m} - \sqrt{m-r} = \sqrt{m} (1 - \sqrt{1 - r/m}) \le \sqrt{m} (r/m + (r/m)^2)/2`.

.. function:: void arb_sqrtpos(arb_t z, const arb_t x, slong prec)

    Sets *z* to the square root of *x*, assuming that *x* represents a
    nonnegative number (i.e. discarding any negative numbers in the input
    interval).

.. function:: void arb_hypot(arb_t z, const arb_t x, const arb_t y, slong prec)

    Sets *z* to `\sqrt{x^2 + y^2}`.

.. function:: void arb_rsqrt(arb_t z, const arb_t x, slong prec)

.. function:: void arb_rsqrt_ui(arb_t z, ulong x, slong prec)

    Sets *z* to the reciprocal square root of *x*, rounded to *prec* bits.
    At high precision, this is faster than computing a square root.

.. function:: void arb_sqrt1pm1(arb_t z, const arb_t x, slong prec)

    Sets `z = \sqrt{1+x}-1`, computed accurately when `x \approx 0`.

.. function:: void arb_root_ui(arb_t z, const arb_t x, ulong k, slong prec)

    Sets *z* to the *k*-th root of *x*, rounded to *prec* bits.
    This function selects between different algorithms. For large *k*,
    it evaluates `\exp(\log(x)/k)`. For small *k*, it uses :func:`arf_root`
    at the midpoint and computes a propagated error bound as follows:
    if input interval is `[m-r, m+r]` with `r \le m`, the error is largest at
    `m-r` where it satisfies

    .. math ::

        m^{1/k} - (m-r)^{1/k} = m^{1/k} [1 - (1-r/m)^{1/k}]

        = m^{1/k} [1 - \exp(\log(1-r/m)/k)]

        \le m^{1/k} \min(1, -\log(1-r/m)/k)

        = m^{1/k} \min(1, \log(1+r/(m-r))/k).

    This is evaluated using :func:`mag_log1p`.

.. function:: void arb_root(arb_t z, const arb_t x, ulong k, slong prec)

    Alias for :func:`arb_root_ui`, provided for backwards compatibility.

.. function:: void arb_sqr(arb_t res, const arb_t x, slong prec)

    Sets *res* to the square of `x`.

.. function:: void _arb_pow_mpz_binexp(arb_t y, const arb_t b, mpz_srcptr e, slong prec)

.. function:: void arb_pow_ui_binexp(arb_t y, const arb_t b, ulong e, slong prec)

    Sets `y = b^e` via binary exponentiation with an initial division if
    `e < 0`. Provided that `b` and `e` are small enough and the exponent is
    positive, the exact power can be computed by setting the precision to
    *ARF_PREC_EXACT*.

    Note that these functions can get slow if the exponent is
    extremely large (in such cases :func:`arb_pow` may be superior).

.. function:: void arb_pow_si(arb_t y, const arb_t b, slong e, slong prec);

    Sets `y = b^e` by wrapping *arb_pow_ui*, doing an initial division if
    `e < 0`.

.. function:: void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, slong prec)

    Sets `y = b^e` by wrapping *arb_pow_si* and *_arb_pow_mpz_binexp*.

.. function:: void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, slong prec)

    Wrapper for *arb_pow_fmpz_binexp*.

.. function:: void arb_pow_ui(arb_t y, const arb_t b, ulong e, slong prec)

    Wrapper for *arb_pow_ui_binexp*.

.. function:: void arb_ui_pow_ui(arb_t y, ulong b, ulong e, slong prec)

.. function:: void arb_si_pow_ui(arb_t y, slong b, ulong e, slong prec)

    Sets `y = b^e` using binary exponentiation. Provided that *b* and *e*
    are small enough and the exponent is positive, the exact power can be
    computed by setting the precision to *ARF_PREC_EXACT*.

.. function:: void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, slong prec)

    Sets `y = b^e`, computed as `y = (b^{1/q})^p` if the denominator of
    `e = p/q` is small, and generally as `y = \exp(e \log b)`.

    Note that this function can get slow if the exponent is
    extremely large (in such cases :func:`arb_pow` may be superior).

.. function:: void arb_pow(arb_t z, const arb_t x, const arb_t y, slong prec)

    Sets `z = x^y`, computed using binary exponentiation if `y` is
    a small exact integer, as `z = (x^{1/2})^{2y}` if `y` is a small exact
    half-integer, and generally as `z = \exp(y \log x)`.

Exponentials and logarithms
-------------------------------------------------------------------------------

.. function:: void arb_log_ui(arb_t z, ulong x, slong prec)

.. function:: void arb_log_fmpz(arb_t z, const fmpz_t x, slong prec)

.. function:: void arb_log_arf(arb_t z, const arf_t x, slong prec)

.. function:: void arb_log(arb_t z, const arb_t x, slong prec)

    Sets `z = \log(x)`.

    At low to medium precision (up to about 4096 bits), :func:`arb_log_arf`
    uses table-based argument reduction and fast Taylor series evaluation
    via :func:`_arb_atan_taylor_rs`. At high precision, it falls back to MPFR.
    The function :func:`arb_log` simply calls :func:`arb_log_arf` with
    the midpoint as input, and separately adds the propagated error.

.. function:: void arb_log_ui_from_prev(arb_t log_k1, ulong k1, arb_t log_k0, ulong k0, slong prec)

    Computes `\log(k_1)`, given `\log(k_0)` where `k_0 < k_1`.
    At high precision, this function uses the formula
    `\log(k_1) = \log(k_0) + 2 \operatorname{atanh}((k_1-k_0)/(k_1+k_0))`,
    evaluating the inverse hyperbolic tangent using binary splitting
    (for best efficiency, `k_0` should be large and `k_1 - k_0` should
    be small). Otherwise, it ignores `\log(k_0)` and evaluates the logarithm
    the usual way.

.. function:: void arb_log1p(arb_t z, const arb_t x, slong prec)

    Sets `z = \log(1+x)`, computed accurately when `x \approx 0`.

.. function:: void arb_log_base_ui(arb_t res, const arb_t x, ulong b, slong prec)

    Sets *res* to `\log_b(x)`. The result is computed exactly when possible.

.. function:: void arb_log_hypot(arb_t res, const arb_t x, const arb_t y, slong prec)

    Sets *res* to `\log(\sqrt{x^2+y^2})`.

.. function:: void arb_exp(arb_t z, const arb_t x, slong prec)

    Sets `z = \exp(x)`. Error propagation is done using the following rule:
    assuming `x = m \pm r`, the error is largest at `m + r`, and we have
    `\exp(m+r) - \exp(m) = \exp(m) (\exp(r)-1) \le r \exp(m+r)`.

.. function:: void arb_expm1(arb_t z, const arb_t x, slong prec)

    Sets `z = \exp(x)-1`, using a more accurate method when `x \approx 0`.

.. function:: void arb_exp_invexp(arb_t z, arb_t w, const arb_t x, slong prec)

    Sets `z = \exp(x)` and `w = \exp(-x)`. The second exponential is computed
    from the first using a division, but propagated error bounds are
    computed separately.

Trigonometric functions
-------------------------------------------------------------------------------

.. function:: void arb_sin(arb_t s, const arb_t x, slong prec)

.. function:: void arb_cos(arb_t c, const arb_t x, slong prec)

.. function:: void arb_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec)

    Sets `s = \sin(x)`, `c = \cos(x)`.

.. function:: void arb_sin_pi(arb_t s, const arb_t x, slong prec)

.. function:: void arb_cos_pi(arb_t c, const arb_t x, slong prec)

.. function:: void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, slong prec)

    Sets `s = \sin(\pi x)`, `c = \cos(\pi x)`.

.. function:: void arb_tan(arb_t y, const arb_t x, slong prec)

    Sets `y = \tan(x) = \sin(x) / \cos(y)`.

.. function:: void arb_cot(arb_t y, const arb_t x, slong prec)

    Sets `y = \cot(x) = \cos(x) / \sin(y)`.

.. function:: void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, slong prec)

.. function:: void arb_sin_pi_fmpq(arb_t s, const fmpq_t x, slong prec)

.. function:: void arb_cos_pi_fmpq(arb_t c, const fmpq_t x, slong prec)

    Sets `s = \sin(\pi x)`, `c = \cos(\pi x)` where `x` is a rational
    number (whose numerator and denominator are assumed to be reduced).
    We first use trigonometric symmetries to reduce the argument to the
    octant `[0, 1/4]`. Then we either multiply by a numerical approximation
    of `\pi` and evaluate the trigonometric function the usual way,
    or we use algebraic methods, depending on which is estimated to be faster.
    Since the argument has been reduced to the first octant, the
    first of these two methods gives full accuracy even if the original
    argument is close to some root other the origin.

.. function:: void arb_tan_pi(arb_t y, const arb_t x, slong prec)

    Sets `y = \tan(\pi x)`.

.. function:: void arb_cot_pi(arb_t y, const arb_t x, slong prec)

    Sets `y = \cot(\pi x)`.

.. function:: void arb_sec(arb_t res, const arb_t x, slong prec)

    Computes `\sec(x) = 1 / \cos(x)`.

.. function:: void arb_csc(arb_t res, const arb_t x, slong prec)

    Computes `\csc(x) = 1 / \sin(x)`.

.. function:: void arb_csc_pi(arb_t res, const arb_t x, slong prec)

    Computes `\csc(\pi x) = 1 / \sin(\pi x)`.

.. function:: void arb_sinc(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{sinc}(x) = \sin(x) / x`.

.. function:: void arb_sinc_pi(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{sinc}(\pi x) = \sin(\pi x) / (\pi x)`.

Inverse trigonometric functions
-------------------------------------------------------------------------------

.. function:: void arb_atan_arf(arb_t z, const arf_t x, slong prec)

.. function:: void arb_atan(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{atan}(x)`.

    At low to medium precision (up to about 4096 bits), :func:`arb_atan_arf`
    uses table-based argument reduction and fast Taylor series evaluation
    via :func:`_arb_atan_taylor_rs`. At high precision, it falls back to MPFR.
    The function :func:`arb_atan` simply calls :func:`arb_atan_arf` with
    the midpoint as input, and separately adds the propagated error.

    The function :func:`arb_atan_arf` uses lookup tables if
    possible, and otherwise falls back to :func:`arb_atan_arf_bb`.

.. function:: void arb_atan2(arb_t z, const arb_t b, const arb_t a, slong prec)

    Sets *r* to an the argument (phase) of the complex number
    `a + bi`, with the branch cut discontinuity on `(-\infty,0]`.
    We define `\operatorname{atan2}(0,0) = 0`, and for `a < 0`,
    `\operatorname{atan2}(0,a) = \pi`.

.. function:: void arb_asin(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{asin}(x) = \operatorname{atan}(x / \sqrt{1-x^2})`.
    If `x` is not contained in the domain `[-1,1]`, the result is an
    indeterminate interval.

.. function:: void arb_acos(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{acos}(x) = \pi/2 - \operatorname{asin}(x)`.
    If `x` is not contained in the domain `[-1,1]`, the result is an
    indeterminate interval.

Hyperbolic functions
-------------------------------------------------------------------------------

.. function:: void arb_sinh(arb_t s, const arb_t x, slong prec)

.. function:: void arb_cosh(arb_t c, const arb_t x, slong prec)

.. function:: void arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec)

    Sets `s = \sinh(x)`, `c = \cosh(x)`. If the midpoint of `x` is close
    to zero and the hyperbolic sine is to be computed,
    evaluates `(e^{2x}\pm1) / (2e^x)` via :func:`arb_expm1`
    to avoid loss of accuracy. Otherwise evaluates `(e^x \pm e^{-x}) / 2`.

.. function:: void arb_tanh(arb_t y, const arb_t x, slong prec)

    Sets `y = \tanh(x) = \sinh(x) / \cosh(x)`, evaluated
    via :func:`arb_expm1` as `\tanh(x) = (e^{2x} - 1) / (e^{2x} + 1)`
    if `|x|` is small, and as
    `\tanh(\pm x) = 1 - 2 e^{\mp 2x} / (1 + e^{\mp 2x})`
    if `|x|` is large.

.. function:: void arb_coth(arb_t y, const arb_t x, slong prec)

    Sets `y = \coth(x) = \cosh(x) / \sinh(x)`, evaluated using
    the same strategy as :func:`arb_tanh`.

.. function:: void arb_sech(arb_t res, const arb_t x, slong prec)

    Computes `\operatorname{sech}(x) = 1 / \cosh(x)`.

.. function:: void arb_csch(arb_t res, const arb_t x, slong prec)

    Computes `\operatorname{csch}(x) = 1 / \sinh(x)`.

Inverse hyperbolic functions
-------------------------------------------------------------------------------

.. function:: void arb_atanh(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{atanh}(x)`.

.. function:: void arb_asinh(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{asinh}(x)`.

.. function:: void arb_acosh(arb_t z, const arb_t x, slong prec)

    Sets `z = \operatorname{acosh}(x)`.
    If `x < 1`, the result is an indeterminate interval.


Constants
-------------------------------------------------------------------------------

The following functions cache the computed values to speed up repeated
calls at the same or lower precision.
For further implementation details, see :ref:`algorithms_constants`.

.. function:: void arb_const_pi(arb_t z, slong prec)

    Computes `\pi`.

.. function:: void arb_const_sqrt_pi(arb_t z, slong prec)

    Computes `\sqrt{\pi}`.

.. function:: void arb_const_log_sqrt2pi(arb_t z, slong prec)

    Computes `\log \sqrt{2 \pi}`.

.. function:: void arb_const_log2(arb_t z, slong prec)

    Computes `\log(2)`.

.. function:: void arb_const_log10(arb_t z, slong prec)

    Computes `\log(10)`.

.. function:: void arb_const_euler(arb_t z, slong prec)

    Computes Euler's constant `\gamma = \lim_{k \rightarrow \infty} (H_k - \log k)`
    where `H_k = 1 + 1/2 + \ldots + 1/k`.

.. function:: void arb_const_catalan(arb_t z, slong prec)

    Computes Catalan's constant `C = \sum_{n=0}^{\infty} (-1)^n / (2n+1)^2`.

.. function:: void arb_const_e(arb_t z, slong prec)

    Computes `e = \exp(1)`.

.. function:: void arb_const_khinchin(arb_t z, slong prec)

    Computes Khinchin's constant `K_0`.

.. function:: void arb_const_glaisher(arb_t z, slong prec)

    Computes the Glaisher-Kinkelin constant `A = \exp(1/12 - \zeta'(-1))`.

.. function:: void arb_const_apery(arb_t z, slong prec)

    Computes Apery's constant `\zeta(3)`.

Lambert W function
-------------------------------------------------------------------------------

.. function:: void arb_lambertw(arb_t res, const arb_t x, int flags, slong prec)

    Computes the Lambert W function, which solves the equation `w e^w = x`.

    The Lambert W function has infinitely many complex branches `W_k(x)`,
    two of which are real on a part of the real line.
    The principal branch `W_0(x)` is selected by setting *flags* to 0, and the
    `W_{-1}` branch is selected by setting *flags* to 1.
    The principal branch is real-valued for `x \ge -1/e`
    (taking values in `[-1,+\infty)`) and the `W_{-1}` branch is real-valued
    for `-1/e \le x < 0` and takes values in `(-\infty,-1]`.
    Elsewhere, the Lambert W function is complex and :func:`acb_lambertw`
    should be used.

    The implementation first computes a floating-point approximation
    heuristically and then computes a rigorously certified enclosure around
    this approximation. Some asymptotic cases are handled specially.
    The algorithm used to compute the Lambert W function is described
    in [Joh2017b]_, which follows the main ideas in [CGHJK1996]_.

Gamma function and factorials
-------------------------------------------------------------------------------

.. function:: void arb_rising_ui(arb_t z, const arb_t x, ulong n, slong prec)
              void arb_rising(arb_t z, const arb_t x, const arb_t n, slong prec)

    Computes the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)`.
    These functions are aliases for :func:`arb_hypgeom_rising_ui`
    and :func:`arb_hypgeom_rising`.

.. function:: void arb_rising_fmpq_ui(arb_t z, const fmpq_t x, ulong n, slong prec)

    Computes the rising factorial `z = x (x+1) (x+2) \cdots (x+n-1)` using
    binary splitting. If the denominator or numerator of *x* is large
    compared to *prec*, it is more efficient to convert *x* to an approximation
    and use :func:`arb_rising_ui`.

.. function :: void arb_rising2_ui(arb_t u, arb_t v, const arb_t x, ulong n, slong prec)

    Letting `u(x) = x (x+1) (x+2) \cdots (x+n-1)`, simultaneously compute
    `u(x)` and `v(x) = u'(x)`.
    This function is a wrapper of :func:`arb_hypgeom_rising_ui_jet`.

.. function:: void arb_fac_ui(arb_t z, ulong n, slong prec)

    Computes the factorial `z = n!` via the gamma function.

.. function:: void arb_doublefac_ui(arb_t z, ulong n, slong prec)

    Computes the double factorial `z = n!!` via the gamma function.

.. function:: void arb_bin_ui(arb_t z, const arb_t n, ulong k, slong prec)

.. function:: void arb_bin_uiui(arb_t z, ulong n, ulong k, slong prec)

    Computes the binomial coefficient `z = {n \choose k}`, via the
    rising factorial as `{n \choose k} = (n-k+1)_k / k!`.

.. function:: void arb_gamma(arb_t z, const arb_t x, slong prec)
              void arb_gamma_fmpq(arb_t z, const fmpq_t x, slong prec)
              void arb_gamma_fmpz(arb_t z, const fmpz_t x, slong prec)

    Computes the gamma function `z = \Gamma(x)`.

    These functions are aliases for :func:`arb_hypgeom_gamma`,
    :func:`arb_hypgeom_gamma_fmpq`, :func:`arb_hypgeom_gamma_fmpz`.

.. function:: void arb_lgamma(arb_t z, const arb_t x, slong prec)

    Computes the logarithmic gamma function `z = \log \Gamma(x)`.
    The complex branch structure is assumed, so if `x \le 0`, the
    result is an indeterminate interval.
    This function is an alias for :func:`arb_hypgeom_lgamma`.

.. function:: void arb_rgamma(arb_t z, const arb_t x, slong prec)

    Computes the reciprocal gamma function `z = 1/\Gamma(x)`,
    avoiding division by zero at the poles of the gamma function.
    This function is an alias for :func:`arb_hypgeom_rgamma`.

.. function:: void arb_digamma(arb_t y, const arb_t x, slong prec)

    Computes the digamma function `z = \psi(x) = (\log \Gamma(x))' = \Gamma'(x) / \Gamma(x)`.


Zeta function
-------------------------------------------------------------------------------

.. function:: void arb_zeta_ui_vec_borwein(arb_ptr z, ulong start, slong num, ulong step, slong prec)

    Evaluates `\zeta(s)` at `\mathrm{num}` consecutive integers *s* beginning
    with *start* and proceeding in increments of *step*.
    Uses Borwein's formula ([Bor2000]_, [GS2003]_),
    implemented to support fast multi-evaluation
    (but also works well for a single *s*).

    Requires `\mathrm{start} \ge 2`. For efficiency, the largest *s*
    should be at most about as
    large as *prec*. Arguments approaching *LONG_MAX* will cause
    overflows.
    One should therefore only use this function for *s* up to about *prec*, and
    then switch to the Euler product.

    The algorithm for single *s* is basically identical to the one used in MPFR
    (see [MPFR2012]_ for a detailed description).
    In particular, we evaluate the sum backwards to avoid storing more than one
    `d_k` coefficient, and use integer arithmetic throughout since it
    is convenient and the terms turn out to be slightly larger than
    `2^\mathrm{prec}`.
    The only numerical error in the main loop comes from the division by `k^s`,
    which adds less than 1 unit of error per term.
    For fast multi-evaluation, we repeatedly divide by `k^{\mathrm{step}}`.
    Each division reduces the input error and adds at most 1 unit of
    additional rounding error, so by induction, the error per term
    is always smaller than 2 units.

.. function:: void arb_zeta_ui_asymp(arb_t x, ulong s, slong prec)

.. function:: void arb_zeta_ui_euler_product(arb_t z, ulong s, slong prec)

    Computes `\zeta(s)` using the Euler product. This is fast only if *s*
    is large compared to the precision. Both methods are trivial wrappers
    for :func:`_acb_dirichlet_euler_product_real_ui`.

.. function:: void arb_zeta_ui_bernoulli(arb_t x, ulong s, slong prec)

    Computes `\zeta(s)` for even *s* via the corresponding Bernoulli number.

.. function:: void arb_zeta_ui_borwein_bsplit(arb_t x, ulong s, slong prec)

    Computes `\zeta(s)` for arbitrary `s \ge 2` using a binary splitting
    implementation of Borwein's algorithm. This has quasilinear complexity
    with respect to the precision (assuming that `s` is fixed).

.. function:: void arb_zeta_ui_vec(arb_ptr x, ulong start, slong num, slong prec)

.. function:: void arb_zeta_ui_vec_even(arb_ptr x, ulong start, slong num, slong prec)

.. function:: void arb_zeta_ui_vec_odd(arb_ptr x, ulong start, slong num, slong prec)

    Computes `\zeta(s)` at *num* consecutive integers (respectively *num*
    even or *num* odd integers) beginning with `s = \mathrm{start} \ge 2`,
    automatically choosing an appropriate algorithm.

.. function:: void arb_zeta_ui(arb_t x, ulong s, slong prec)

    Computes `\zeta(s)` for nonnegative integer `s \ne 1`, automatically
    choosing an appropriate algorithm. This function is
    intended for numerical evaluation of isolated zeta values; for
    multi-evaluation, the vector versions are more efficient.

.. function:: void arb_zeta(arb_t z, const arb_t s, slong prec)

    Sets *z* to the value of the Riemann zeta function `\zeta(s)`.

    For computing derivatives with respect to `s`,
    use :func:`arb_poly_zeta_series`.

.. function:: void arb_hurwitz_zeta(arb_t z, const arb_t s, const arb_t a, slong prec)

    Sets *z* to the value of the Hurwitz zeta function `\zeta(s,a)`.

    For computing derivatives with respect to `s`,
    use :func:`arb_poly_zeta_series`.

Bernoulli numbers and polynomials
-------------------------------------------------------------------------------

.. function:: void arb_bernoulli_ui(arb_t b, ulong n, slong prec)

.. function:: void arb_bernoulli_fmpz(arb_t b, const fmpz_t n, slong prec)

    Sets `b` to the numerical value of the Bernoulli number `B_n`
    approximated to *prec* bits.

    The internal precision is increased automatically to give an accurate
    result. Note that, with huge *fmpz* input, the output will have a huge
    exponent and evaluation will accordingly be slower.

    A single division from the exact fraction of `B_n` is used if this value
    is in the global cache or the exact numerator roughly is larger than
    *prec* bits. Otherwise, the Riemann zeta function is used
    (see :func:`arb_bernoulli_ui_zeta`).

    This function reads `B_n` from the global cache
    if the number is already cached, but does not automatically extend
    the cache by itself.

.. function:: void arb_bernoulli_ui_zeta(arb_t b, ulong n, slong prec)

    Sets `b` to the numerical value of `B_n` accurate to *prec* bits,
    computed using the formula `B_{2n} = (-1)^{n+1} 2 (2n)! \zeta(2n) / (2 \pi)^n`.

    To avoid potential infinite recursion, we explicitly call the
    Euler product implementation of the zeta function.
    This method will only give high accuracy if the precision is small
    enough compared to `n` for the Euler product to converge rapidly.

.. function:: void arb_bernoulli_poly_ui(arb_t res, ulong n, const arb_t x, slong prec)

    Sets *res* to the value of the Bernoulli polynomial `B_n(x)`.

    Warning: this function is only fast if either *n* or *x* is a small integer.

    This function reads Bernoulli numbers from the global cache if they
    are already cached, but does not automatically extend the cache by itself.

.. function:: void arb_power_sum_vec(arb_ptr res, const arb_t a, const arb_t b, slong len, slong prec)

    For *n* from 0 to *len* - 1, sets entry *n* in the output vector *res* to

    .. math ::

        S_n(a,b) = \frac{1}{n+1}\left(B_{n+1}(b) - B_{n+1}(a)\right)

    where `B_n(x)` is a Bernoulli polynomial. If *a* and *b* are integers
    and `b \ge a`, this is equivalent to

    .. math ::

        S_n(a,b) = \sum_{k=a}^{b-1} k^n.

    The computation uses the generating function for Bernoulli polynomials.

Polylogarithms
-------------------------------------------------------------------------------

.. function:: void arb_polylog(arb_t w, const arb_t s, const arb_t z, slong prec)

.. function:: void arb_polylog_si(arb_t w, slong s, const arb_t z, slong prec)

    Sets *w* to the polylogarithm `\operatorname{Li}_s(z)`.

Other special functions
-------------------------------------------------------------------------------

.. function:: void arb_fib_fmpz(arb_t z, const fmpz_t n, slong prec)

.. function:: void arb_fib_ui(arb_t z, ulong n, slong prec)

    Computes the Fibonacci number `F_n`. Uses the binary squaring
    algorithm described in [Tak2000]_.
    Provided that *n* is small enough, an exact Fibonacci number can be
    computed by setting the precision to *ARF_PREC_EXACT*.

.. function:: void arb_agm(arb_t z, const arb_t x, const arb_t y, slong prec)

    Sets *z* to the arithmetic-geometric mean of *x* and *y*.

.. function:: void arb_chebyshev_t_ui(arb_t a, ulong n, const arb_t x, slong prec)

.. function:: void arb_chebyshev_u_ui(arb_t a, ulong n, const arb_t x, slong prec)

    Evaluates the Chebyshev polynomial of the first kind `a = T_n(x)`
    or the Chebyshev polynomial of the second kind `a = U_n(x)`.

.. function:: void arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)

.. function:: void arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)

    Simultaneously evaluates `a = T_n(x), b = T_{n-1}(x)` or
    `a = U_n(x), b = U_{n-1}(x)`.
    Aliasing between *a*, *b* and *x* is not permitted.

.. function:: void arb_bell_sum_bsplit(arb_t res, const fmpz_t n, const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec)

.. function:: void arb_bell_sum_taylor(arb_t res, const fmpz_t n, const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec)

    Helper functions for Bell numbers, evaluating the sum
    `\sum_{k=a}^{b-1} k^n / k!`. If *mmag* is non-NULL, it may be used
    to indicate that the target error tolerance should be
    `2^{mmag - prec}`.

.. function:: void arb_bell_fmpz(arb_t res, const fmpz_t n, slong prec)

.. function:: void arb_bell_ui(arb_t res, ulong n, slong prec)

    Sets *res* to the Bell number `B_n`. If the number is too large to
    fit exactly in *prec* bits, a numerical approximation is computed
    efficiently.

    The algorithm to compute Bell numbers, including error analysis,
    is described in detail in [Joh2015]_.

.. function:: void arb_euler_number_fmpz(arb_t res, const fmpz_t n, slong prec)
              void arb_euler_number_ui(arb_t res, ulong n, slong prec)

    Sets *res* to the Euler number `E_n`, which is defined by
    the exponential generating function `1 / \cosh(x)`.
    The result will be exact if `E_n` is exactly representable
    at the requested precision.

.. function:: void arb_fmpz_euler_number_ui_multi_mod(fmpz_t res, ulong n, double alpha)
              void arb_fmpz_euler_number_ui(fmpz_t res, ulong n)

    Computes the Euler number `E_n` as an exact integer. The default
    algorithm uses a table lookup, the Dirichlet beta function or a
    hybrid modular algorithm depending on the size of *n*.
    The *multi_mod* algorithm accepts a tuning parameter *alpha* which
    can be set to a negative value to use defaults.

.. function:: void arb_partitions_fmpz(arb_t res, const fmpz_t n, slong prec)

.. function:: void arb_partitions_ui(arb_t res, ulong n, slong prec)

    Sets *res* to the partition function `p(n)`.
    When *n* is large and `\log_2 p(n)` is more than twice *prec*,
    the leading term in the Hardy-Ramanujan asymptotic series is used
    together with an error bound. Otherwise, the exact value is computed
    and rounded.

.. function:: void arb_primorial_nth_ui(arb_t res, ulong n, slong prec)

    Sets *res* to the *nth* primorial, defined as the product of the
    first *n* prime numbers. The running time is quasilinear in *n*.

.. function:: void arb_primorial_ui(arb_t res, ulong n, slong prec)

    Sets *res* to the primorial defined as the product of the positive
    integers up to and including *n*. The running time is quasilinear in *n*.

Internals for computing elementary functions
-------------------------------------------------------------------------------

.. function:: void _arb_atan_taylor_naive(mp_ptr y, mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N, int alternating)

.. function:: void _arb_atan_taylor_rs(mp_ptr y, mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N, int alternating)

    Computes an approximation of `y = \sum_{k=0}^{N-1} x^{2k+1} / (2k+1)`
    (if *alternating* is 0) or `y = \sum_{k=0}^{N-1} (-1)^k x^{2k+1} / (2k+1)`
    (if *alternating* is 1). Used internally for computing arctangents
    and logarithms. The *naive* version uses the forward recurrence, and the
    *rs* version uses a division-avoiding rectangular splitting scheme.

    Requires `N \le 255`, `0 \le x \le 1/16`, and *xn* positive.
    The input *x* and output *y* are fixed-point numbers with *xn* fractional
    limbs. A bound for the ulp error is written to *error*.

.. function:: void _arb_exp_taylor_naive(mp_ptr y, mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N)

.. function:: void _arb_exp_taylor_rs(mp_ptr y, mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N)

    Computes an approximation of `y = \sum_{k=0}^{N-1} x^k / k!`. Used internally
    for computing exponentials. The *naive* version uses the forward recurrence,
    and the *rs* version uses a division-avoiding rectangular splitting scheme.

    Requires `N \le 287`, `0 \le x \le 1/16`, and *xn* positive.
    The input *x* is a fixed-point number with *xn* fractional
    limbs, and the output *y* is a fixed-point number with *xn* fractional
    limbs plus one extra limb for the integer part of the result.

    A bound for the ulp error is written to *error*.

.. function:: void _arb_sin_cos_taylor_naive(mp_ptr ysin, mp_ptr ycos, mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N)

.. function:: void _arb_sin_cos_taylor_rs(mp_ptr ysin, mp_ptr ycos, mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N, int sinonly, int alternating)

    Computes approximations of `y_s = \sum_{k=0}^{N-1} (-1)^k x^{2k+1} / (2k+1)!`
    and `y_c = \sum_{k=0}^{N-1} (-1)^k x^{2k} / (2k)!`.
    Used internally for computing sines and cosines. The *naive* version uses
    the forward recurrence, and the *rs* version uses a division-avoiding
    rectangular splitting scheme.

    Requires `N \le 143`, `0 \le x \le 1/16`, and *xn* positive.
    The input *x* and outputs *ysin*, *ycos* are fixed-point numbers with
    *xn* fractional limbs. A bound for the ulp error is written to *error*.

    If *sinonly* is 1, only the sine is computed; if *sinonly* is 0
    both the sine and cosine are computed.
    To compute sin and cos, *alternating* should be 1. If *alternating* is 0,
    the hyperbolic sine is computed (this is currently only intended to
    be used together with *sinonly*).

.. function:: int _arb_get_mpn_fixed_mod_log2(mp_ptr w, fmpz_t q, mp_limb_t * error, const arf_t x, mp_size_t wn)

    Attempts to write `w = x - q \log(2)` with `0 \le w < \log(2)`, where *w*
    is a fixed-point number with *wn* limbs and ulp error *error*.
    Returns success.

.. function:: int _arb_get_mpn_fixed_mod_pi4(mp_ptr w, fmpz_t q, int * octant, mp_limb_t * error, const arf_t x, mp_size_t wn)

    Attempts to write `w = |x| - q \pi/4` with `0 \le w < \pi/4`, where *w*
    is a fixed-point number with *wn* limbs and ulp error *error*.
    Returns success.

    The value of *q* mod 8 is written to *octant*. The output variable *q*
    can be NULL, in which case the full value of *q* is not stored.

.. function:: slong _arb_exp_taylor_bound(slong mag, slong prec)

    Returns *n* such that
    `\left|\sum_{k=n}^{\infty} x^k / k!\right| \le 2^{-\mathrm{prec}}`,
    assuming `|x| \le 2^{\mathrm{mag}} \le 1/4`.

.. function:: void arb_exp_arf_bb(arb_t z, const arf_t x, slong prec, int m1)

    Computes the exponential function using the bit-burst algorithm.
    If *m1* is nonzero, the exponential function minus one is computed
    accurately.

    Aborts if *x* is extremely small or large (where another algorithm
    should be used).

    For large *x*, repeated halving is used. In fact, we always
    do argument reduction until `|x|` is smaller than about `2^{-d}`
    where `d \approx 16` to speed up convergence. If `|x| \approx 2^m`,
    we thus need about `m+d` squarings.

    Computing `\log(2)` costs roughly 100-200 multiplications, so is not
    usually worth the effort at very high precision. However, this function
    could be improved by using `\log(2)` based reduction at precision low
    enough that the value can be assumed to be cached.

.. function:: void _arb_exp_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)

.. function:: void _arb_exp_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)

    Computes *T*, *Q* and *Qexp* such that
    `T / (Q 2^{\text{Qexp}}) = \sum_{k=1}^N (x/2^r)^k/k!` using binary splitting.
    Note that the sum is taken to *N* inclusive and omits the constant term.

    The *powtab* version precomputes a table of powers of *x*,
    resulting in slightly higher memory usage but better speed. For best
    efficiency, *N* should have many trailing zero bits.

.. function:: void arb_exp_arf_rs_generic(arb_t res, const arf_t x, slong prec, int minus_one)

    Computes the exponential function using a generic version of the rectangular
    splitting strategy, intended for intermediate precision.

.. function:: void _arb_atan_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)

.. function:: void _arb_atan_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp, const fmpz_t x, flint_bitcnt_t r, slong N)

    Computes *T*, *Q* and *Qexp* such that
    `T / (Q 2^{\text{Qexp}}) = \sum_{k=1}^N (-1)^k (x/2^r)^{2k} / (2k+1)`
    using binary splitting.
    Note that the sum is taken to *N* inclusive, omits the linear term,
    and requires a final multiplication by `(x/2^r)` to give the
    true series for atan.

    The *powtab* version precomputes a table of powers of *x*,
    resulting in slightly higher memory usage but better speed. For best
    efficiency, *N* should have many trailing zero bits.

.. function:: void arb_atan_arf_bb(arb_t z, const arf_t x, slong prec)

    Computes the arctangent of *x*.
    Initially, the argument-halving formula

    .. math ::

        \operatorname{atan}(x) = 2 \operatorname{atan}\left(\frac{x}{1+\sqrt{1+x^2}}\right)

    is applied up to 8 times to get a small argument.
    Then a version of the bit-burst algorithm is used.
    The functional equation

    .. math ::

        \operatorname{atan}(x) = \operatorname{atan}(p/q) +
            \operatorname{atan}(w),
            \quad w = \frac{qx-p}{px+q},
            \quad p = \lfloor qx \rfloor

    is applied repeatedly instead of integrating a differential
    equation for the arctangent, as this appears to be more efficient.

.. function:: void arb_sin_cos_arf_generic(arb_t s, arb_t c, const arf_t x, slong prec)

    Computes the sine and cosine of *x* using a generic strategy.
    This function gets called internally by the main sin and cos functions
    when the precision for argument reduction or series evaluation
    based on lookup tables is exhausted.

    This function first performs a cheap test to see if
    `|x| < \pi / 2 - \varepsilon`. If the test fails, it uses
    `\pi` to reduce the argument to the first octant,
    and then evaluates the sin and cos functions recursively (this call cannot
    result in infinite recursion).

    If no argument reduction is needed, this function uses a generic version
    of the rectangular splitting algorithm if the precision is not too high,
    and otherwise invokes the asymptotically fast bit-burst algorithm.

.. function:: void arb_sin_cos_arf_bb(arb_t s, arb_t c, const arf_t x, slong prec)

    Computes the sine and cosine of *x* using the bit-burst algorithm.
    It is required that `|x| < \pi / 2` (this is not checked).

.. function:: void arb_sin_cos_wide(arb_t s, arb_t c, const arb_t x, slong prec)

    Computes an accurate enclosure (with both endpoints optimal to within
    about `2^{-30}` as afforded by the radius format) of the range of
    sine and cosine on a given wide interval. The computation is done
    by evaluating the sine and cosine at the interval endpoints and
    determining whether peaks of -1 or 1 occur between the endpoints.
    The interval is then converted back to a ball.

    The internal computations are done with doubles, using a simple
    floating-point algorithm to approximate the sine and cosine. It is
    easy to see that the cumulative errors in this algorithm add up to
    less than `2^{-30}`, with the dominant source of error being a single
    approximate reduction by `\pi/2`. This reduction is done safely using
    doubles up to a magnitude of about `2^{20}`. For larger arguments, a
    slower reduction using :type:`arb_t` arithmetic is done as a
    preprocessing step.

.. function:: void arb_sin_cos_generic(arb_t s, arb_t c, const arb_t x, slong prec)

    Computes the sine and cosine of *x* by taking care of various special
    cases and computing the propagated error before calling
    :func:`arb_sin_cos_arf_generic`. This is used as a fallback inside
    :func:`arb_sin_cos` to take care of all cases without a fast
    path in that function.

Vector functions
-------------------------------------------------------------------------------

.. function:: void _arb_vec_zero(arb_ptr vec, slong n)

    Sets all entries in *vec* to zero.

.. function:: int _arb_vec_is_zero(arb_srcptr vec, slong len)

    Returns nonzero iff all entries in *x* are zero.

.. function:: int _arb_vec_is_finite(arb_srcptr x, slong len)

    Returns nonzero iff all entries in *x* certainly are finite.

.. function:: void _arb_vec_set(arb_ptr res, arb_srcptr vec, slong len)

    Sets *res* to a copy of *vec*.

.. function:: void _arb_vec_set_round(arb_ptr res, arb_srcptr vec, slong len, slong prec)

    Sets *res* to a copy of *vec*, rounding each entry to *prec* bits.

.. function:: void _arb_vec_swap(arb_ptr vec1, arb_ptr vec2, slong len)

    Swaps the entries of *vec1* and *vec2*.

.. function:: void _arb_vec_neg(arb_ptr B, arb_srcptr A, slong n)

.. function:: void _arb_vec_sub(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)

.. function:: void _arb_vec_add(arb_ptr C, arb_srcptr A, arb_srcptr B, slong n, slong prec)

.. function:: void _arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)

.. function:: void _arb_vec_scalar_div(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)

.. function:: void _arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec, slong len, const fmpz_t c, slong prec)

.. function:: void _arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c)

.. function:: void _arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec, slong len, const arb_t c, slong prec)

   Performs the respective scalar operation elementwise.

.. function:: void _arb_vec_get_mag(mag_t bound, arb_srcptr vec, slong len, slong prec)

    Sets *bound* to an upper bound for the entries in *vec*.

.. function:: slong _arb_vec_bits(arb_srcptr x, slong len)

    Returns the maximum of :func:`arb_bits` for all entries in *vec*.

.. function:: void _arb_vec_set_powers(arb_ptr xs, const arb_t x, slong len, slong prec)

    Sets *xs* to the powers `1, x, x^2, \ldots, x^{len-1}`.

.. function:: void _arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, slong len)

.. function:: void _arb_vec_add_error_mag_vec(arb_ptr res, mag_srcptr err, slong len)

    Adds the magnitude of each entry in *err* to the radius of the
    corresponding entry in *res*.

.. function:: void _arb_vec_indeterminate(arb_ptr vec, slong len)

    Applies :func:`arb_indeterminate` elementwise.

.. function:: void _arb_vec_trim(arb_ptr res, arb_srcptr vec, slong len)

    Applies :func:`arb_trim` elementwise.

.. function:: int _arb_vec_get_unique_fmpz_vec(fmpz * res,  arb_srcptr vec, slong len)

    Calls :func:`arb_get_unique_fmpz` elementwise and returns nonzero if
    all entries can be rounded uniquely to integers. If any entry in *vec*
    cannot be rounded uniquely to an integer, returns zero.

