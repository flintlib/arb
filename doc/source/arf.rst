.. _arf:

**arf.h** -- arbitrary-precision floating-point numbers
===============================================================================

A variable of type :type:`arf_t` holds an arbitrary-precision binary
floating-point number, i.e. a rational number of the form
`x \times 2^y` where `x, y \in \mathbb{Z}` and `x` is odd;
or one of the special values zero, plus infinity, minus infinity,
or NaN (not-a-number).

The *exponent* of a finite and nonzero floating-point number can be
defined in different
ways: for example, as the component *y* above, or as the unique
integer *e* such that
`x \times 2^y = m \times 2^e` where `1/2 \le |m| < 1`.
The internal representation of an :type:`arf_t` stores the
exponent in the latter format.

The conventions for special values largely follow those of the
IEEE floating-point standard. At the moment, there is no support
for negative zero, unsigned infinity, or a NaN with a payload, though
some these might be added in the future.

Except where otherwise noted, the output of an operation is the
floating-point number obtained by taking the inputs as exact numbers,
in principle carrying out the operation exactly, and rounding the
resulting real number to the nearest representable floating-point
number whose mantissa has at most the specified number of bits, in
the specified direction of rounding. Some operations are always
or optionally done exactly.

The :type:`arf_t` type is almost identical semantically to
the legacy :type:`fmpr_t` type, but uses a more efficient
internal representation.
The most significant differences that the user
has to be aware of are:

* The mantissa is no longer represented as a FLINT :type:`fmpz`, and the
  internal exponent points to the top of the binary expansion of the mantissa
  instead of of the bottom. Code designed to manipulate components of an
  :type:`fmpr_t` directly can be ported to the :type:`arf_t` type 
  by making use of :func:`arf_get_fmpz_2exp` and :func:`arf_set_fmpz_2exp`.

* Some :type:`arf_t` functions return an :type:`int`
  indicating whether a result is inexact, whereas the corresponding
  :type:`fmpr_t` functions return a :type:`long` encoding the relative
  exponent of the error.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: arf_struct

.. type:: arf_t

    An :type:`arf_struct` contains four words: an :type:`fmpz` exponent (*exp*),
    a *size* field tracking the number of limbs used (one bit of this
    field is also used for the sign of the number), and two more words.
    The last two words hold the value directly if there are at most two limbs,
    and otherwise contain one *alloc* field (tracking the total number of
    allocated limbs, not all of which might be used) and a pointer to
    the actual limbs.
    Thus, up to 128 bits on a 64-bit machine and 64 bits on a 32-bit machine,
    no space outside of the :type:`arf_struct` is used.

    An :type:`arf_t` is defined as an array of length one of type
    :type:`arf_struct`, permitting an :type:`arf_t` to be passed by reference.

.. type:: arf_rnd_t

    Specifies the rounding mode for the result of an approximate operation.

.. macro:: ARF_RND_DOWN

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards zero.

.. macro:: ARF_RND_UP

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction away from zero.

.. macro:: ARF_RND_FLOOR

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards minus infinity.

.. macro:: ARF_RND_CEIL

    Specifies that the result of an operation should be rounded to the
    nearest representable number in the direction towards plus infinity.

.. macro:: ARF_RND_NEAR

    Specifies that the result of an operation should be rounded to the
    nearest representable number, rounding to an odd mantissa if there is a tie
    between two values. *Warning*: this rounding mode is currently
    not implemented (except for a few conversions functions where this 
    stated explicitly).

.. macro:: ARF_PREC_EXACT

    If passed as the precision parameter to a function, indicates that no
    rounding is to be performed. This must only be used when it is known
    that the result of the operation can be represented exactly and fits
    in memory (the typical use case is working with small integer values).
    Note that, for example, adding two numbers whose exponents are far
    apart can easily produce an exact result that is far too large to
    store in memory.

Memory management
-------------------------------------------------------------------------------

.. function:: void arf_init(arf_t x)

    Initializes the variable *x* for use. Its value is set to zero.

.. function:: void arf_clear(arf_t x)

    Clears the variable *x*, freeing or recycling its allocated memory.

Special values
-------------------------------------------------------------------------------

.. function:: void arf_zero(arf_t x)

.. function:: void arf_one(arf_t x)

.. function:: void arf_pos_inf(arf_t x)

.. function:: void arf_neg_inf(arf_t x)

.. function:: void arf_nan(arf_t x)

    Sets *x* respectively to 0, 1, `+\infty`, `-\infty`, NaN.

.. function:: int arf_is_zero(const arf_t x)

.. function:: int arf_is_one(const arf_t x)

.. function:: int arf_is_pos_inf(const arf_t x)

.. function:: int arf_is_neg_inf(const arf_t x)

.. function:: int arf_is_nan(const arf_t x)

    Returns nonzero iff *x* respectively equals 0, 1, `+\infty`, `-\infty`, NaN.

.. function:: int arf_is_inf(const arf_t x)

    Returns nonzero iff *x* equals either `+\infty` or `-\infty`.

.. function:: int arf_is_normal(const arf_t x)

    Returns nonzero iff *x* is a finite, nonzero floating-point value, i.e.
    not one of the special values 0, `+\infty`, `-\infty`, NaN.

.. function:: int arf_is_special(const arf_t x)

    Returns nonzero iff *x* is one of the special values
    0, `+\infty`, `-\infty`, NaN, i.e. not a finite, nonzero
    floating-point value.

.. function:: int arf_is_finite(arf_t x)

    Returns nonzero iff *x* is a finite floating-point value,
    i.e. not one of the values `+\infty`, `-\infty`, NaN.
    (Note that this is not equivalent to the negation of
    :func:`arf_is_inf`.)


Assignment, rounding and conversions
-------------------------------------------------------------------------------

.. function:: void arf_set(arf_t y, const arf_t x)

.. function:: void arf_set_mpz(arf_t y, const mpz_t x)

.. function:: void arf_set_fmpz(arf_t y, const fmpz_t x)

.. function:: void arf_set_ui(arf_t y, ulong x)

.. function:: void arf_set_si(arf_t y, long x)

.. function:: void arf_set_mpfr(arf_t y, const mpfr_t x)

.. function:: void arf_set_fmpr(arf_t y, const fmpr_t x)

.. function:: void arf_set_d(arf_t y, double x)

    Sets *y* exactly to *x*.

.. function:: void arf_swap(arf_t y, arf_t x)

    Swaps *y* and *x* efficiently.

.. function:: void arf_init_set_ui(arf_t y, ulong x)

.. function:: void arf_init_set_si(arf_t y, long x)

    Initialises *y* and sets it to *x* in a single operation.

.. function:: int arf_set_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)

.. function:: int arf_set_round_si(arf_t x, long v, long prec, arf_rnd_t rnd)

.. function:: int arf_set_round_ui(arf_t x, ulong v, long prec, arf_rnd_t rnd)

.. function:: int arf_set_round_mpz(arf_t y, const mpz_t x, long prec, arf_rnd_t rnd)

.. function:: int arf_set_round_fmpz(arf_t y, const fmpz_t x, long prec, arf_rnd_t rnd)

    Sets *y* to *x*, rounded to *prec* bits in the direction
    specified by *rnd*.

.. function:: void arf_set_si_2exp_si(arf_t y, long m, long e)

.. function:: void arf_set_ui_2exp_si(arf_t y, ulong m, long e)

.. function:: void arf_set_fmpz_2exp(arf_t y, const fmpz_t m, const fmpz_t e)

    Sets *y* to `m \times 2^e`.

.. function:: int arf_set_round_fmpz_2exp(arf_t y, const fmpz_t x, const fmpz_t e, long prec, arf_rnd_t rnd)

    Sets *y* to `x \times 2^e`, rounded to *prec* bits in the direction
    specified by *rnd*.

.. function:: void arf_get_fmpz_2exp(fmpz_t m, fmpz_t e, const arf_t x)

    Sets *m* and *e* to the unique integers such that
    `x = m \times 2^e` and *m* is odd,
    provided that *x* is a nonzero finite fraction.
    If *x* is zero, both *m* and *e* are set to zero. If *x* is
    infinite or NaN, the result is undefined.

.. function:: double arf_get_d(const arf_t x, arf_rnd_t rnd)

    Returns *x* rounded to a double in the direction specified by *rnd*.

.. function:: void arf_get_fmpr(fmpr_t y, const arf_t x)

    Sets *y* exactly to *x*.

.. function:: int arf_get_mpfr(mpfr_t y, const arf_t x, mpfr_rnd_t rnd)

    Sets the MPFR variable *y* to the value of *x*. If the precision of *x*
    is too small to allow *y* to be represented exactly, it is rounded in
    the specified MPFR rounding mode. The return value (-1, 0 or 1)
    indicates the direction of rounding, following the convention
    of the MPFR library.

.. function:: void arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)

    Sets *z* to *x* rounded to the nearest integer in the direction
    specified by *rnd*. If rnd is *ARF_RND_NEAR*, rounds to the nearest
    even integer in case of a tie. Aborts if *x* is infinite, NaN or if the
    exponent is unreasonably large.

.. function:: long arf_get_si(const arf_t x, arf_rnd_t rnd)

    Returns *x* rounded to the nearest integer in the direction specified by
    *rnd*. If *rnd* is *ARF_RND_NEAR*, rounds to the nearest even integer
    in case of a tie. Aborts if *x* is infinite, NaN, or the value is
    too large to fit in a long.

.. function:: int arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e)

.. function:: int arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, long e)

    Converts *x* to a mantissa with predetermined exponent, i.e. computes
    an integer *y* such that `y \times 2^e \approx x`, truncating if necessary.
    Returns 0 if exact and 1 if truncation occurred.

.. function:: void arf_floor(arf_t y, const arf_t x)

.. function:: void arf_ceil(arf_t y, const arf_t x)

    Sets *y* to `\lfloor x \rfloor` and `\lceil x \rceil` respectively.
    The result is always represented exactly, requiring no more bits to
    store than the input. To round the result to a floating-point number
    with a lower precision, call :func:`arf_set_round` afterwards.

Comparisons and bounds
-------------------------------------------------------------------------------

.. function:: int arf_equal(const arf_t x, const arf_t y)

.. function:: int arf_equal_si(const arf_t x, long y)

    Returns nonzero iff *x* and *y* are exactly equal. This function does
    not treat NaN specially, i.e. NaN compares as equal to itself.

.. function:: int arf_cmp(const arf_t x, const arf_t y)

    Returns negative, zero, or positive, depending on whether *x* is
    respectively smaller, equal, or greater compared to *y*.
    Comparison with NaN is undefined.

.. function:: int arf_cmpabs(const arf_t x, const arf_t y)

.. function:: int arf_cmpabs_ui(const arf_t x, ulong y)

.. function:: int arf_cmpabs_mag(const arf_t x, const mag_t y)

    Compares the absolute values of *x* and *y*.

.. function:: int arf_cmp_2exp_si(const arf_t x, long e)

.. function:: int arf_cmpabs_2exp_si(const arf_t x, long e)

    Compares *x* (respectively its absolute value) with `2^e`.

.. function:: int arf_sgn(const arf_t x)

    Returns `-1`, `0` or `+1` according to the sign of *x*. The sign
    of NaN is undefined.

.. function:: void arf_min(arf_t z, const arf_t a, const arf_t b)

.. function:: void arf_max(arf_t z, const arf_t a, const arf_t b)

    Sets *z* respectively to the minimum and the maximum of *a* and *b*.

.. function:: long arf_bits(const arf_t x)

    Returns the number of bits needed to represent the absolute value
    of the mantissa of *x*, i.e. the minimum precision sufficient to represent
    *x* exactly. Returns 0 if *x* is a special value.

.. function:: int arf_is_int(const arf_t x)

    Returns nonzero iff *x* is integer-valued.

.. function:: int arf_is_int_2exp_si(const arf_t x, long e)

    Returns nonzero iff *x* equals `n 2^e` for some integer *n*.

.. function:: void arf_abs_bound_lt_2exp_fmpz(fmpz_t b, const arf_t x)

    Sets *b* to the smallest integer such that `|x| < 2^b`.
    If *x* is zero, infinity or NaN, the result is undefined.

.. function:: void arf_abs_bound_le_2exp_fmpz(fmpz_t b, const arf_t x)

    Sets *b* to the smallest integer such that `|x| \le 2^b`.
    If *x* is zero, infinity or NaN, the result is undefined.

.. function:: long arf_abs_bound_lt_2exp_si(const arf_t x)

    Returns the smallest integer *b* such that `|x| < 2^b`, clamping
    the result to lie between -*ARF_PREC_EXACT* and *ARF_PREC_EXACT*
    inclusive. If *x* is zero, -*ARF_PREC_EXACT* is returned,
    and if *x* is infinity or NaN, *ARF_PREC_EXACT* is returned.

Magnitude functions
-------------------------------------------------------------------------------

.. function:: void arf_get_mag(mag_t y, const arf_t x)

    Sets *y* to an upper bound for the absolute value of *x*.

.. function:: void arf_get_mag_lower(mag_t y, const arf_t x)

    Sets *y* to a lower bound for the absolute value of *x*.

.. function:: void arf_set_mag(arf_t y, const mag_t x)

    Sets *y* to *x*.

.. function:: void mag_init_set_arf(mag_t y, const arf_t x)

    Initializes *y* and sets it to an upper bound for *x*.

.. function:: void mag_fast_init_set_arf(mag_t y, const arf_t x)

    Initializes *y* and sets it to an upper bound for *x*.
    Assumes that the exponent of *y* is small.

.. function:: void arf_mag_set_ulp(mag_t z, const arf_t y, long prec)

    Sets *z* to the magnitude of the unit in the last place (ulp) of *y*
    at precision *prec*.

.. function:: void arf_mag_add_ulp(mag_t z, const mag_t x, const arf_t y, long prec)

    Sets *z* to an upper bound for the sum of *x* and the
    magnitude of the unit in the last place (ulp) of *y*
    at precision *prec*.

.. function:: void arf_mag_fast_add_ulp(mag_t z, const mag_t x, const arf_t y, long prec)

    Sets *z* to an upper bound for the sum of *x* and the
    magnitude of the unit in the last place (ulp) of *y*
    at precision *prec*. Assumes that all exponents are small.

Shallow assignment
-------------------------------------------------------------------------------

.. function:: void arf_init_set_shallow(arf_t z, const arf_t x)

.. function:: void arf_init_set_mag_shallow(arf_t z, const mag_t x)

    Initializes *z* to a shallow copy of *x*. A shallow copy just involves
    copying struct data (no heap allocation is performed).

    The target variable *z* may not be cleared or modified in any way (it can
    only be used as constant input to functions), and may not be used after
    *x* has been cleared. Moreover, after *x* has been assigned shallowly
    to *z*, no modification of *x* is permitted as long as *z* is in use.

.. function:: void arf_init_neg_shallow(arf_t z, const arf_t x)

.. function:: void arf_init_neg_mag_shallow(arf_t z, const mag_t x)

    Initializes *z* shallowly to the negation of *x*.

Random number generation
-------------------------------------------------------------------------------

.. function:: void arf_randtest(arf_t x, flint_rand_t state, long bits, long mag_bits)

    Generates a finite random number whose mantissa has precision at most
    *bits* and whose exponent has at most *mag_bits* bits. The
    values are distributed non-uniformly: special bit patterns are generated
    with high probability in order to allow the test code to exercise corner
    cases.

.. function:: void arf_randtest_not_zero(arf_t x, flint_rand_t state, long bits, long mag_bits)

    Identical to :func:`arf_randtest`, except that zero is never produced
    as an output.

.. function:: void arf_randtest_special(arf_t x, flint_rand_t state, long bits, long mag_bits)

    Indentical to :func:`arf_randtest`, except that the output occasionally
    is set to an infinity or NaN.

Input and output
-------------------------------------------------------------------------------

.. function:: void arf_debug(const arf_t x)

    Prints information about the internal representation of *x*.

.. function:: void arf_print(const arf_t x)

    Prints *x* as an integer mantissa and exponent.

.. function:: void arf_printd(const arf_t y, long d)

    Prints *x* as a decimal floating-point number, rounding to *d* digits.
    This function is currently implemented using MPFR,
    and does not support large exponents.

Addition and multiplication
-------------------------------------------------------------------------------

.. function:: void arf_abs(arf_t y, const arf_t x)

    Sets *y* to the absolute value of *x*.

.. function:: void arf_neg(arf_t y, const arf_t x)

    Sets `y = -x` exactly.

.. function:: int arf_neg_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)

    Sets `y = -x`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

.. function:: void arf_mul_2exp_si(arf_t y, const arf_t x, long e)

.. function:: void arf_mul_2exp_fmpz(arf_t y, const arf_t x, const fmpz_t e)

    Sets `y = x 2^e` exactly.

.. function:: int arf_mul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_mul_ui(arf_t z, const arf_t x, ulong y, long prec, arf_rnd_t rnd)

.. function:: int arf_mul_si(arf_t z, const arf_t x, long y, long prec, arf_rnd_t rnd)

.. function:: int arf_mul_mpz(arf_t z, const arf_t x, const mpz_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_mul_fmpz(arf_t z, const arf_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

    Sets `z = x \times y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

.. function:: int arf_add(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_add_si(arf_t z, const arf_t x, long y, long prec, arf_rnd_t rnd)

.. function:: int arf_add_ui(arf_t z, const arf_t x, ulong y, long prec, arf_rnd_t rnd)

.. function:: int arf_add_fmpz(arf_t z, const arf_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

    Sets `z = x + y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

.. function:: int arf_add_fmpz_2exp(arf_t z, const arf_t x, const fmpz_t y, const fmpz_t e, long prec, arf_rnd_t rnd)

    Sets `z = x + y 2^e`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

.. function:: int arf_sub(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_sub_si(arf_t z, const arf_t x, long y, long prec, arf_rnd_t rnd)

.. function:: int arf_sub_ui(arf_t z, const arf_t x, ulong y, long prec, arf_rnd_t rnd)

.. function:: int arf_sub_fmpz(arf_t z, const arf_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

    Sets `z = x - y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

.. function:: int arf_addmul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_addmul_ui(arf_t z, const arf_t x, ulong y, long prec, arf_rnd_t rnd)

.. function:: int arf_addmul_si(arf_t z, const arf_t x, long y, long prec, arf_rnd_t rnd)

.. function:: int arf_addmul_mpz(arf_t z, const arf_t x, const mpz_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_addmul_fmpz(arf_t z, const arf_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

    Sets `z = z + x \times y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

.. function:: int arf_submul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_submul_ui(arf_t z, const arf_t x, ulong y, long prec, arf_rnd_t rnd)

.. function:: int arf_submul_si(arf_t z, const arf_t x, long y, long prec, arf_rnd_t rnd)

.. function:: int arf_submul_mpz(arf_t z, const arf_t x, const mpz_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_submul_fmpz(arf_t z, const arf_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

    Sets `z = z - x \times y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact.

Summation
-------------------------------------------------------------------------------

.. function:: int arf_sum(arf_t s, arf_srcptr terms, long len, long prec, arf_rnd_t rnd)

    Sets *s* to the sum of the array *terms* of length *len*, rounded to
    *prec* bits in the direction specified by *rnd*. The sum is computed as if
    done without any intermediate rounding error, with only a single rounding
    applied to the final result. Unlike repeated calls to :func:`arf_add` with
    infinite precision, this function does not overflow if the magnitudes of
    the terms are far apart. Warning: this function is implemented naively,
    and the running time is quadratic with respect to *len* in the worst case.

Division
-------------------------------------------------------------------------------

.. function:: int arf_div(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_div_ui(arf_t z, const arf_t x, ulong y, long prec, arf_rnd_t rnd)

.. function:: int arf_ui_div(arf_t z, ulong x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_div_si(arf_t z, const arf_t x, long y, long prec, arf_rnd_t rnd)

.. function:: int arf_si_div(arf_t z, long x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_div_fmpz(arf_t z, const arf_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_fmpz_div(arf_t z, const fmpz_t x, const arf_t y, long prec, arf_rnd_t rnd)

.. function:: int arf_fmpz_div_fmpz(arf_t z, const fmpz_t x, const fmpz_t y, long prec, arf_rnd_t rnd)

    Sets `z = x / y`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact. The result is NaN if *y* is zero.

Square roots
-------------------------------------------------------------------------------

.. function:: int arf_sqrt(arf_t z, const arf_t x, long prec, arf_rnd_t rnd)

.. function:: int arf_sqrt_ui(arf_t z, ulong x, long prec, arf_rnd_t rnd)

.. function:: int arf_sqrt_fmpz(arf_t z, const fmpz_t x, long prec, arf_rnd_t rnd)

    Sets `z = \sqrt{x}`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact. The result is NaN if *x* is negative.

.. function:: int arf_rsqrt(arf_t z, const arf_t x, long prec, arf_rnd_t rnd)

    Sets `z = 1/\sqrt{x}`, rounded to *prec* bits in the direction specified by *rnd*,
    returning nonzero iff the operation is inexact. The result is NaN if *x* is
    negative, and `+\infty` if *x* is zero.

Complex arithmetic
-------------------------------------------------------------------------------

.. function:: int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, long prec, arf_rnd_t rnd)

.. function:: int arf_complex_mul_fallback(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, long prec, arf_rnd_t rnd)

    Computes the complex product `e + fi = (a + bi)(c + di)`, rounding both
    `e` and `f` correctly to *prec* bits in the direction specified by *rnd*.
    The first bit in the return code indicates inexactness of `e`, and the
    second bit indicates inexactness of `f`.

    If any of the components *a*, *b*, *c*, *d* is zero, two real
    multiplications and no additions are done. This convention is used even
    if any other part contains an infinity or NaN, and the behavior
    with infinite/NaN input is defined accordingly.

    The *fallback* version is implemented naively, for testing purposes.
    No squaring optimization is implemented.

.. function:: int arf_complex_sqr(arf_t e, arf_t f, const arf_t a, const arf_t b, long prec, arf_rnd_t rnd)

    Computes the complex square `e + fi = (a + bi)^2`. This function has
    identical semantics to :func:`arf_complex_mul` (with `c = a, b = d`),
    but is faster.

