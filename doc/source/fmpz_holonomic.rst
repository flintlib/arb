fmpz_holonomic.h -- holonomic sequences and functions
===============================================================================

This module provides support for holonomic sequences and functions,
which are represented via annihilating operators
`A = \sum_{i=0}^r a_i D^i` where `a_i \in \mathbb{Z}[x]`.

Depending on context, `D` is one of the following:

A forward shift operator `D f(x) = f(x+1)`, making `A` a difference operator
which defines a family of sequences `c(x)` satisfying the difference equation

.. math ::

    a_0(x) c(x) + a_1(x) c(x+1) + \ldots + a_r(x) c(x+r) = 0.

A differential operator `D f(x) = f'(x)`, making `A` a differential operator
which defines a family of functions `f(x)` satisfying the differential
equation

.. math ::

    a_0(x) f(x) + a_1(x) f'(x) + \ldots + a_r(x) f^{(r)}(x) = 0.

A specific sequence of function is defined by fixing `r - 1` initial values
(derivatives), where `r` is called the *order* of the operator `A`.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_holonomic_struct

.. type:: fmpz_holonomic_t

    Represents a holonomic difference or differential operator as
    a sequence of coefficients of type *fmpz_poly_struct*.

    An *fmpz_holonomic_t* is defined as an array of length one of type
    *fmpz_holonomic_t*, permitting an *fmpz_holonomic_t* to
    be passed by reference.


Memory management
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_init(fmpz_holonomic_t op)

    Initializes *op* for use, setting it to the zero operator (note that
    this is not a valid operator for input to most functions).


.. function:: void fmpz_holonomic_clear(fmpz_holonomic_t op)

    Clears *op*, deallocating all coefficients and the coefficient array.


.. function:: void fmpz_holonomic_fit_length(fmpz_holonomic_t op, long len)

    Makes sure that the coefficient array of the operator contains at
    least *len* initialized polynomials.


.. function:: void _fmpz_holonomic_set_length(fmpz_holonomic_t op, long len)

    Directly changes the length of the operator, without allocating or
    deallocating coefficients. The value shold not exceed the allocation length.


Basic operations
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_set(fmpz_holonomic_t rop, const fmpz_holonomic_t op)

    Sets `op` to a copy of `op`.

.. function:: void fmpz_holonomic_one(fmpz_holonomic_t op)

    Sets *op* to the constant operator 1. As a differential operator, this
    annihilates the zero function. As a difference operator, it annihilates
    the zero sequence.

.. function:: void fmpz_holonomic_randtest(fmpz_holonomic_t op, flint_rand_t state, long r, long d, long b)

    Sets *op* to a random nonzero operator of order at most *r*, degree
    at most *d*, and with coefficients at most *b* bits in size. The operator
    is guaranteed to have a nonzero leading coefficient, but otherwise
    will not be normalised.

Properties
-------------------------------------------------------------------------------

.. function:: long fmpz_holonomic_order(const fmpz_holonomic_t op)

    Returns the order *r* of *op*.

.. function:: long fmpz_holonomic_degree(const fmpz_holonomic_t op)

    Returns the degree *d* of *op*, defined as the highest degree of
    all its coefficients.

.. function:: int fmpz_holonomic_seq_is_constant(const fmpz_holonomic_t op)

    Returns nonzero if *op* is zero-order, i.e. annihilates constant
    sequences.

.. function:: int fmpz_holonomic_seq_is_cfinite(const fmpz_holonomic_t op)

    Returns nonzero if *op* has constant coefficients, i.e. annihilates
    C-finite sequences.

.. function:: int fmpz_holonomic_seq_is_hypgeom(const fmpz_holonomic_t op)

    Return nonzero if *op* is first-order, i.e. annihilates hypergeometric
    sequences.

Input and output
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_print(const fmpz_holonomic_t op, const char * x, const char * d)

    Prints a pretty representation of *op*, using the string *x* for the
    variable of the coefficient polynomials, and using the string *d* for
    the differential or difference operator.


Normalisation
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_normalise_leading(fmpz_holonomic_t op)

    Normalises *op* by removing leading zero coefficients.

.. function:: void fmpz_holonomic_normalise_sign(fmpz_holonomic_t op)

    Normalises *op* by making the leading coefficient of the leading
    polynomial positive.

.. function:: void fmpz_holonomic_normalise_content(fmpz_holonomic_t op)

    Normalises *op* by dividing out the content, i.e. the greatest common
    divisor, of all the coefficients.

.. function:: void fmpz_holonomic_seq_normalise_trailing(fmpz_holonomic_t op)

    Normalises *op* as a difference operator by removing trailing
    zero coefficients. This requires shifting the higher-order coefficients
    to compensate.

.. function:: void fmpz_holonomic_seq_normalise(fmpz_holonomic_t op)

    Normalises *op* as a difference operator by removing leading and trailing
    zero coefficients, removing the content, and making the leading
    polynomial positive.


Shifting
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_shift_fmpz(fmpz_holonomic_t res, const fmpz_holonomic_t op, const fmpz_t s)

.. function:: void fmpz_holonomic_shift_fmpq(fmpz_holonomic_t res, const fmpz_holonomic_t op, const fmpq_t s)

.. function:: void fmpz_holonomic_shift_si(fmpz_holonomic_t res, const fmpz_holonomic_t op, long s)

    Given an operator *op* annihilating a function or sequence `f(x)`,
    sets *res* to an annihilator of the shifted function
    or sequence `f(x+s)`.


Special sequences
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_seq_set_const(fmpz_holonomic_t op)

    Sets *op* to an annihilator of the constant sequence `c, c, c, \ldots`
    where `c` is arbitrary.

.. function:: void fmpz_holonomic_seq_set_fmpz_pow(fmpz_holonomic_t op, const fmpz_t c)

.. function:: void fmpz_holonomic_seq_set_fmpq_pow(fmpz_holonomic_t op, const fmpq_t c)

    Sets *op* to an annihilator of the sequence `c, c^2, c^3, \ldots`.

.. function:: void fmpz_holonomic_seq_set_factorial(fmpz_holonomic_t op)

    Sets *op* to an annihilator of the sequence of factorials `n!`.

.. function:: void fmpz_holonomic_seq_set_harmonic(fmpz_holonomic_t op)

    Sets *op* to an annihilator of the sequence of harmonic numbers
    `H_n = 1 + 1/2 + 1/3 + \ldots + 1/n`.

.. function:: void fmpz_holonomic_seq_set_fibonacci(fmpz_holonomic_t op)

    Sets *op* to an annihilator of the sequence of Fibonacci numbers `F_n`.


Special functions
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_fun_set_pow_fmpq(fmpz_holonomic_t op, const fmpq_t c)

.. function:: void fmpz_holonomic_fun_set_pow_fmpz(fmpz_holonomic_t op, const fmpz_t c)

    Sets *op* to a differential operator annihilating the power function `x^c`.

.. function:: void fmpz_holonomic_fun_set_exp(fmpz_holonomic_t op)

    Sets *op* to a differential operator annihilating the exponential function `e^x`.

.. function:: void fmpz_holonomic_fun_set_sin_cos(fmpz_holonomic_t op)

    Sets *op* to a differential operator annihilating the sine
    and cosine functions `\sin(x)` and `\cos(x)`.

.. function:: void fmpz_holonomic_fun_set_log(fmpz_holonomic_t op)

    Sets *op* to a differential operator annihilating the
    natural logarithm `\log(x)`.

.. function:: void fmpz_holonomic_fun_set_atan(fmpz_holonomic_t op)

    Sets *op* to a differential operator annihilating the
    inverse tangent function `\operatorname{atan}(x)`.

.. function:: void fmpz_holonomic_fun_set_asin_acos(fmpz_holonomic_t op)

    Sets *op* to a differential operator annihilating the inverse
    sine and cosine functions `\operatorname{asin}(x)` and `\operatorname{acos}(x)`.

.. function:: void fmpz_holonomic_fun_set_erf(fmpz_holonomic_t op)

    Sets *op* to a differential operator annihilating the error function `\operatorname{erf}(x)`.


Sequence transformations
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_seq_mul(fmpz_holonomic_t res, const fmpz_holonomic_t op1, const fmpz_holonomic_t op2)

    Given annihilators *op1* and *op2* of sequences `a_0, a_1, \ldots` and
    `b_0, b_1, \ldots`, sets *res* to an annihilator of the sequence
    `a_0 b_0, a_1 b_1, \ldots`.
    This function currently requires *op1* and *op2* to be hypergeometric
    (i.e. of order 1).

.. function:: void fmpz_holonomic_seq_pow_si(fmpz_holonomic_t res, const fmpz_holonomic_t op, long e)

    Given an annihilator *op* of the sequence `c_0, c_1, c_2, \ldots`, outputs
    an annihilator of the sequence `c_0^e, c_1^e, c_2^e, \ldots`.
    This function currently requires *op* to be hypergeometric
    (i.e. of order 1).

.. function:: void fmpz_holonomic_seq_reverse(fmpz_holonomic_t res, const fmpz_holonomic_t op)

    Given an annihilator *op* of the sequence `c(n)`, sets *res* to an
    annihilator of the sequence `c(-n)`.

.. function:: void fmpz_holonomic_seq_section(fmpz_holonomic_t res, const fmpz_holonomic_t op, long m)

    Given an annihilator *op* of the sequence `c(n)` and an integer constant *m*,
    sets *res* to an annihilator of the sequence `c(mn)`.
    The constant *m* can be zero or negative.


Taylor series
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_get_series(fmpz_holonomic_t re, const fmpz_holonomic_t de)

    Given a differential operator *de* annihilating some function `f(x)`,
    sets *re* to a difference operator annihilating the coefficients `c_k`
    in the Taylor series at `x = 0`,

    .. math ::

        f(x) = \sum_{k=0}^{\infty} c_k x^k.


Sequence evaluation
-------------------------------------------------------------------------------

.. function:: void fmpz_holonomic_forward_fmpz_mat(fmpz_mat_t M, fmpz_t Q, const fmpz_holonomic_t op, long start, long n)

    Let *op* be an operator of order *r* annihilating a sequence `c(k)`.
    This function computes an *r* by *r* integer matrix `M` and a
    denominator `Q` such that, for any initial values
    `c(s), c(s+1), \ldots, c(s+r-1)` where `s` is given by *start*,

    .. math ::

        Q \, \begin{pmatrix} c(s+n) \\ c(s+n+1) \\ \vdots \\ c(s+n+r-1) \end{pmatrix}
        = 
        M \, \begin{pmatrix} c(s) \\ c(s+1) \\ \vdots \\ c(s+r-1) \end{pmatrix}.

    The output is computed by multiplying together successive companion
    matrices, using binary splitting to balance the sizes of the subproducts.

    Some special cases are handled more efficiently. In particular,
    if *op* has constant coefficients, matrix exponentiation is used.

    In general, no attempt is made to divide out content from `M` and `Q`.
    If `Q` is zero, the leading coefficient of *op* has a root somewhere
    among the evaluation points, making the sequence undefined from
    that point onward.

.. function:: void fmpz_holonomic_forward_fmprb_mat(fmprb_mat_t M, fmprb_t Q, const fmpz_holonomic_t op, long start, long n, long prec)

    Equivalent to the *fmpz_mat* version, but truncates large entries.

.. function:: void fmpz_holonomic_get_nth_fmpz(fmpz_t res, const fmpz_holonomic_t op, const fmpz * initial, long n0, long n)

.. function:: void fmpz_holonomic_get_nth_fmpq(fmpq_t res, const fmpz_holonomic_t op, const fmpq * initial, long n0, long n)

    Computes element `c(n)` in the sequence annihilated by the
    difference operator *op*, given the
    initial values `c(n_0), c(n_1), \ldots, c(n_0+r-1)` where
    `r` is the order of *op*.
    The *fmpz* version assumes that the result is actually an integer.


