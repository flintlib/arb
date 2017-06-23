.. _arb_fmpz_poly:

**arb_fmpz_poly.h** -- extra methods for integer polynomials
=========================================================================================

This module provides methods for FLINT polynomials
with integer and rational coefficients
(:type:`fmpz_poly_t`) and (:type:`fmpq_poly_t`)
requiring use of Arb real or complex numbers.

Some methods output real or complex numbers while others
use real and complex numbers internally to produce an exact result.
This module also contains some useful helper functions not specifically
related to real and complex numbers that are useful elsewhere in Arb.

Note that methods that combine Arb *polynomials* and FLINT polynomials
are found in the respective Arb polynomial modules, such as
:func:`arb_poly_set_fmpz_poly`
and :func:`arb_poly_get_unique_fmpz_poly`.

Evaluation
-------------------------------------------------------------------------------

.. function:: void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * poly, slong len, const arb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

.. function:: void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * poly, slong len, const acb_t x, slong prec)

.. function:: void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)

    Evaluates *poly* (given by a polynomial object or an array with *len* coefficients)
    at the given real or complex number, respectively using Horner's rule, rectangular
    splitting, or a default algorithm choice.

Special polynomials
-------------------------------------------------------------------------------

Note: see also the methods available in FLINT (e.g. for cyclotomic polynomials).

.. function:: void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n)

    Sets *res* to the minimal polynomial of the Gaussian periods
    `\sum_{a \in H} \zeta^a` where `\zeta = \exp(2 \pi i / q)`
    and *H* are the cosets of the subgroups of order `d = (q - 1) / n` of
    `(\mathbb{Z}/q\mathbb{Z})^{\times}`.
    The resulting polynomial has degree *n*.
    When `d = 1`, the result is the cyclotomic polynomial `\Phi_q`.

    The implementation assumes that *q* is prime, and that *n* is a divisor of
    `q - 1` such that *n* is coprime with *d*. If any condition is not met,
    *res* is set to the zero polynomial.

    This method provides a fast (in practice) way to
    construct finite field extensions of prescribed degree.
    If *q* satisfies the conditions stated above and `(q-1)/f` additionally
    is coprime with *n*, where *f* is the multiplicative order of *p* mod *q*, then
    the Gaussian period minimal polynomial is irreducible over
    `\operatorname{GF}(p)` [CP2005]_.

