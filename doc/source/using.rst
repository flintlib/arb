.. _using:

Using ball arithmetic
===============================================================================

Ball semantics
-------------------------------------------------------------------------------

Let `f : A \to B` be a function.
A ball implementation of `f` is a function `F` that maps subsets of `A`
to subsets of `B` subject to the following rule:

    For all `x \in X \subseteq A`,
    we have `f(x) \in F(X) \subseteq B`.

In other words, `F(X)` is an *enclosure* for the set `\{f(x) : x \in X\}`.
This rule is sometimes called the *inclusion principle*.

Throughout the documentation (except where otherwise noted),
we will simply write `f(x)` instead of `F(X)`
when describing ball implementations of pointwise-defined mathematical
functions, understanding that the input is a set of point values and that
the output is an enclosure.

General subsets of `\mathbb{R}` are not possible to
represent on a computer. Instead, we work with subsets of the form
`[m \pm r] = [m-r, m+r]` where the midpoint *m* and radius *r* are binary
floating-point numbers, i.e. numbers of the form `u 2^v` with `u, v \in \mathbb{Z}`
(to make this scheme complete,
we also need to adjoin the special floating-point
values `-\infty`, `+\infty` and `\operatorname{NaN}`).

Given a ball `[m \pm r]` with `m \in \mathbb{R}` (not necessarily a
floating-point number),
we can always round *m* to a nearby floating-point number that has at most
most *prec* bits in the component *u*,
and add an upper bound for the rounding error to *r*.
In Arb, ball functions that take a *prec* argument as input always
round their output to *prec* bits.
Some functions are always exact, and thus do not take a *prec* argument.

While the internal representation uses binary floating-point numbers,
it is usually preferable to print numbers in decimal. The binary-to-decimal
conversion generally requires rounding, so a printout
`[m \pm r]` will typically have a slightly larger *r*
than the actual radius that is stored (see :func:`arb_get_str` for details).

Predicates
-------------------------------------------------------------------------------

A ball implementation of a predicate 
`f : \mathbb{R} \to \{\operatorname{True}, \operatorname{False}\}`
would need to be able to return a third logical value indicating
that the result could be either True or False.
In most cases, predicates in Arb are implemented as 
functions that return the *int* value 1 to indicate that the
result certainly is True, and the *int* value 0 to indicate
that the result could be either True or False.
To test whether a predicate certainly is False, the user must
test whether the negated predicate certainly is True.

For example, the following code would *not* be correct in general:

.. code-block:: c

    if (arb_is_positive(x))
    {
        ...  /* do things assuming that x > 0 */
    }
    else
    {
        ...  /* do things assuming that x <= 0 */
    }

Instead, the following can be used:

.. code-block:: c

    if (arb_is_positive(x))
    {
        ...  /* do things assuming that x > 0 */
    }
    else if (arb_is_nonpositive(x))
    {
        ...  /* do things assuming that x <= 0 */
    }
    else
    {
        ... /* do things assuming that the sign of x is unknown */
    }

Likewise, we will write `x \le y` in mathematical notation with the meaning
that `x \le y` holds for all `x \in X, y \in Y` where `X` and `Y` are balls.

Note that that some predicates naturally act on the input *viewed as balls*,
rather than being ball implementations of pointwise predicates.
For example :func:`arb_contains_zero` tests whether the input ball
contains the point zero, and this is an exact test as such.

A worked example: the sine function
-------------------------------------------------------------------------------

We implement the function `\sin(x)` naively using
the Taylor series `\sum_{k=0}^{\infty} (-1)^k x^{2k+1} / (2k+1)!`
and :type:`arb_t` arithmetic (see :ref:`arb`).
Since there are infinitely many terms, we need to split the series
in two parts: a finite sum that can be evaluated directly, and
a tail that has to be bounded.

We stop as soon as we reach a term `t` bounded by `|t| \le 2^{-prec} < 1`.
The terms are alternating and must have decreasing magnitude
from that point, so the tail of the series
is bounded by `|t|`. We add this magnitude to the radius
of the output. Since ball arithmetic automatically bounds the numerical errors
resulting from all arithmetic operations, the output *res* is a
ball guaranteed to contain `\sin(x)`.

.. code-block:: c

    #include "arb.h"

    void arb_sin_naive(arb_t res, const arb_t x, slong prec)
    {
        arb_t s, t, u, tol;
        slong k;
        arb_init(s); arb_init(t); arb_init(u); arb_init(tol);

        arb_one(tol);
        arb_mul_2exp_si(tol, tol, -prec);  /* tol = 2^-prec */

        for (k = 0; ; k++)
        {
            arb_pow_ui(t, x, 2 * k + 1, prec);
            arb_fac_ui(u, 2 * k + 1, prec);
            arb_div(t, t, u, prec);  /* t = x^(2k+1) / (2k+1)! */

            arb_abs(u, t);
            if (arb_le(u, tol))   /* if |t| <= 2^-prec */
            {
                arb_add_error(s, u);    /* add |t| to the radius and stop */
                break;
            }

            if (k % 2 == 0)
                arb_add(s, s, t, prec);
            else
                arb_sub(s, s, t, prec);

        }

        arb_set(res, s);
        arb_clear(s); arb_clear(t); arb_clear(u); arb_clear(tol);
    }

This algorithm is naive, because the Taylor series is slow to converge
and suffers from catastrophic cancellation when `|x|` is large
(we could also improve the efficiency of the code slightly by
computing the terms using recurrence relations instead of
computing `x^k` and `k!` from scratch each iteration).

As a test, we compute `\sin(2016.1)`.
The largest term in the Taylor series for `\sin(x)` reaches
a magnitude of about `x^x / x!`, or about `10^{873}` in this case.
Therefore, we need over 873 digits (about 3000 bits) of precision
to overcome the catastrophic cancellation and determine
the result with sufficient accuracy to tell whether it is positive
or negative.

.. code-block:: c

    int main()
    {
        arb_t x, y;
        slong prec;
        arb_init(x); arb_init(y);

        for (prec = 64; ; prec *= 2)
        {
            arb_set_str(x, "2016.1", prec);
            arb_sin_naive(y, x, prec);
            printf("Using %5ld bits, sin(x) = ", prec);
            arb_printn(y, 10, 0); printf("\n");
            if (!arb_contains_zero(y))  /* stopping condition */
                break;
        }

        arb_clear(x); arb_clear(y);
    }

The program produces the following output::

    Using    64 bits, sin(x) = [+/- 2.67e+859]
    Using   128 bits, sin(x) = [+/- 1.30e+840]
    Using   256 bits, sin(x) = [+/- 3.60e+801]
    Using   512 bits, sin(x) = [+/- 3.01e+724]
    Using  1024 bits, sin(x) = [+/- 2.18e+570]
    Using  2048 bits, sin(x) = [+/- 1.22e+262]
    Using  4096 bits, sin(x) = [-0.7190842207 +/- 1.20e-11]

As an exercise, the reader may improve the naive algorithm by making it
subtract a well-chosen multiple of `2 \pi` from `x` before invoking
the Taylor series (hint: use :func:`arb_const_pi`, :func:`arb_div`
and :func:`arf_get_fmpz`).
If done correctly, 64 bits of precision should be more than enough to
compute `\sin(2016.1)`, and with minor adjustments
to the code, the user should be able to compute
`\sin(\exp(2016.1))` quite easily as well.

This example illustrates how ball arithmetic can be used to perform
nontrivial calculations. To evaluate an infinite series, the user
needs to know how to bound the tail of the series, but everything
else is automatic.
When evaluating a finite formula that can be expressed
completely using built-in functions, all error bounding is automatic
from the point of view of the user.
In particular, the :func:`arb_sin` method should be used to compute the sine
of a real number; it uses a much more efficient algorithm
than the naive code above.

This example also illustrates the "guess-and-verify" paradigm:
instead of determining *a priori* the floating-point precision necessary
to get a correct result, we *guess* some initial precision, use ball arithmetic
to *verify* that the result is accurate enough, and restart with
higher precision (or signal failure) if it is not.

If we think of rounding errors as essentially random processes,
then a floating-point computation is analogous to a
*Monte Carlo algorithm*. Using ball arithmetic to get a verified result
effectively turns it into the analog of a *Las Vegas algorithm*,
which is a randomized algorithm that always gives a correct result if it terminates, but
may fail to terminate (alternatively, instead of actually looping forever,
it might signal failure after a certain number of iterations).

The loop will fail to terminate if we attempt to determine the sign of
`\sin(\pi)`::

    Using    64 bits, sin(x) = [+/- 3.96e-18]
    Using   128 bits, sin(x) = [+/- 2.17e-37]
    Using   256 bits, sin(x) = [+/- 6.10e-76]
    Using   512 bits, sin(x) = [+/- 5.13e-153]
    Using  1024 bits, sin(x) = [+/- 4.01e-307]
    Using  2048 bits, sin(x) = [+/- 2.13e-615]
    Using  4096 bits, sin(x) = [+/- 6.85e-1232]
    Using  8192 bits, sin(x) = [+/- 6.46e-2465]
    Using 16384 bits, sin(x) = [+/- 5.09e-4931]
    Using 32768 bits, sin(x) = [+/- 5.41e-9863]
    ...

The sign of a nonzero real number can be
decided by computing it to sufficiently high accuracy, but the sign
of an expression that is exactly equal to zero cannot be decided
by a numerical computation unless the entire computation
happens to be exact (in this example, we could use the :func:`arb_sin_pi` 
function which computes `\sin(\pi x)` in one step, with the input `x = 1`).

It is up to the user to implement a stopping criterion appropriate for
the circumstances of a given application. For example, breaking
when it is clear that `|\sin(x)| < 10^{-10000}` would allow the program
to terminate and convey some meaningful information about the input `x = \pi`,
though this would not constitute a mathematical proof that
`\sin(\pi) = 0`.

More on precision and accuracy
-------------------------------------------------------------------------------

Basic operations on balls such addition and multiplication only
involve a single floating-point operation on the midpoint.
The effect of the *prec* argument is then obvious.
More complicated functions are computed by performing a long sequence
of arithmetic operations, each of which requires a rounding
and also propagates the error accumulated from previous operations.

The *prec* argument essentially controls the internal working precision
for each step. A higher higher or lower precision might be used internally
in order to try to achieve an accuracy of *prec* bits.
To complicate things further, many algorithms require
approximation steps (such as truncation of infinite series)
that depend on the precision in a more subtle way.
As a result, the relation between *prec* and the accuracy of
the output is not always easy to predict.

(To be expanded.)
