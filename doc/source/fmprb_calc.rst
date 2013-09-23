.. _fmprb-calc:

**fmprb_calc.h** -- calculus with real-valued functions
===============================================================================

This module provides functions for operations of calculus
over the real numbers (intended to include root-finding,
optimization, integration, and so on). It is planned that the module
will include two types of algorithms:

* Interval algorithms that give provably correct results. An example
  would be numerical integration on an interval by dividing the
  interval into small balls and evaluating the function
  on each ball, giving rigorous upper and lower bounds.
* Conventional numerical algorithms that use heuristics
  to estimate the accuracy of a result, without guaranteeing
  that it is correct. An example would be numerical integration
  based on pointwise evaluation, where the error is estimated
  by comparing the results with two different sets of evaluation
  points. Ball arithmetic then still tracks the accuracy
  of the function evaluations.

Any algorithms of the second kind will be clearly
marked as such.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmprb_calc_func_t

    Typedef for a pointer to a function with signature

    .. code ::

        int func(fmprb_ptr out, const fmprb_t inp, void * param, long order, long prec)

    implementing a univariate real function `f(x)`.
    When called, *func* should write to *out* the first *order*
    coefficients in the Taylor series expansion of `f(x)` at the point *inp*,
    evaluated at a precision of *prec* bits.
    The *param* argument may be used to pass through
    additional parameters to the function.
    The return value is reserved for future use as an
    error code. It can be assumed that *out* and *inp* are not
    aliased and that *order* is positive.

.. macro:: FMPRB_CALC_SUCCESS

    Return value indicating that an operation is successful.

.. macro:: FMPRB_CALC_IMPRECISE_INPUT

    Return value indicating that the input to a function probably needs
    to be computed more accurately.

.. macro:: FMPRB_CALC_NO_CONVERGENCE

    Return value indicating that an algorithm has failed to convergence,
    possibly due to the problem not having a solution, the algorithm
    not being applicable, or the precision being insufficient

Debugging
-------------------------------------------------------------------------------

.. var:: int fmprb_calc_verbose

    If set, enables printing information about the calculation
    to standard output.


Root-finding
-------------------------------------------------------------------------------

.. function:: void fmprb_calc_newton_conv_factor(fmpr_t conv_factor, fmprb_calc_func_t func, void * param, const fmprb_t conv_region, long prec)

    Given an interval `I` specified by *conv_region*, evaluates a bound
    for `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    where `f` is the function specified by *func* and *param*.
    The bound is obtained by evaluating `f'(I)` and `f''(I)` directly.
    If `f` is ill-conditioned, `I` may need to be extremely precise in
    order to get an effective, finite bound for *C*.

.. function:: int fmprb_calc_newton_step(fmprb_t xnew, fmprb_calc_func_t func, void * param, const fmprb_t x, const fmprb_t conv_region, const fmpr_t conv_factor, long prec)

    Performs a single step with an interval version of Newton's method.
    The input consists of the function `f` specified
    by *func* and *param*, a ball `x = [m-r, m+r]` known
    to contain a single root of `f`, a ball `I` (*conv_region*)
    containing `x` with an associated bound (*conv_factor*) for
    `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    and a working precision *prec*.

    The Newton update consists of setting
    `x' = [m'-r', m'+r']` where `m' = m - f(m) / f'(m)`
    and `r' = C r^2`. The expression `m - f(m) / f'(m)` is evaluated
    using ball arithmetic at a working precision of *prec* bits, and the
    rounding error during this evaluation is accounted for in the output.
    We now check that `x' \in I` and `r' < r`. If both conditions are
    satisfied, we set *xnew* to `x'` and return *FMPRB_CALC_SUCCESS*.
    If either condition fails, we set *xnew* to `x` and return
    *FMPRB_CALC_NO_CONVERGENCE*, indicating that no progress
    is made.

.. function:: int fmprb_calc_refine_root_newton(fmprb_t r, fmprb_calc_func_t func, void * param, const fmprb_t start, const fmprb_t conv_region, const fmpr_t conv_factor, long eval_extra_prec, long prec)

    Refines a precise estimate of a single root of a function
    to high precision by performing several Newton steps, using
    nearly optimally chosen doubling precision steps.

    The inputs are defined as for *fmprb_calc_newton_step*, except for
    the precision parameters: *prec* is the target accuracy and
    *eval_extra_prec* is the estimated number of guard bits that need
    to be added to evaluate the function accurately close to the root
    (for example, if the function is a polynomial with large coefficients
    of alternating signs and Horner's rule is used to evaluate it,
    the extra precision should typically be approximately
    the bit size of the coefficients).

    This function returns *FMPRB_CALC_SUCCESS* if all attempted
    Newton steps are successful (note that this does not guarantee
    that the computed root is accurate to *prec* bits, which has
    to be verified by the user), only that it is more accurate
    than the starting ball.

    On failure, *FMPRB_CALC_IMPRECISE_INPUT*
    or *FMPRB_CALC_NO_CONVERGENCE* may be returned. In this case, *r*
    is set to a ball for the root which is valid but likely
    does have full accuracy (it can possibly just be equal
    to the starting ball).

