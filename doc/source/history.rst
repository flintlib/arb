History and changes
===============================================================================

For more details, view the detailed commit log
in the git repository https://github.com/fredrik-johansson/arb

* 2013-01-26 - version 0.4

  * much faster fmpr_mul, fmprb_mul and set_round, resulting in general speed improvements
  * code for computing the complex Hurwitz zeta function with derivatives
  * fixed and documented error bounds for hypergeometric series
  * better algorithm for series evaluation of the gamma function at a rational point
  * much faster generation of Bernoulli numbers
  * complex log, exp, pow, trigonometric functions (currently based on MPFR)
  * complex nth roots via Newton iteration
  * added code for arithmetic on fmpcb_polys
  * code for computing Khinchin's constant
  * code for rising factorials of polynomials or power series
  * faster sin_cos
  * better div_2expm1
  * many other new helper functions
  * improved thread safety
  * more test code for core operations

* 2012-11-07 - version 0.3

  * converted documentation to sphinx
  * new module fmpcb for ball interval arithmetic over the complex numbers

    * conversions, utility functions and arithmetic operations

  * new module fmpcb_mat for matrices over the complex numbers

    * conversions, utility functions and arithmetic operations
    * multiplication, LU decomposition, solving, inverse and determinant

  * new module fmpcb_poly for polynomials over the complex numbers

    * root isolation for complex polynomials

  * new module fmpz_holonomic for functions/sequences
    defined by linear differential/difference equations
    with polynomial coefficients

    * functions for creating various special sequences and functions
    * some closure properties for sequences
    * Taylor series expansion for differential equations
    * computing the nth entry of a sequence using binary splitting
    * computing the nth entry mod p using fast multipoint evaluation

  * generic binary splitting code with automatic error bounding is now
    used for evaluating hypergeometric series
  * matrix powering
  * various other helper functions

* 2012-09-29 - version 0.2

  * code for computing the gamma function (Karatsuba, Stirling's series)
  * rising factorials
  * fast exp_series using Newton iteration
  * improved multiplication of small polynomials by using classical multiplication
  * implemented error propagation for square roots
  * polynomial division (Newton-based)
  * polynomial evaluation (Horner) and composition (divide-and-conquer)
  * product trees, fast multipoint evaluation and interpolation (various algorithms)
  * power series composition (Horner, Brent-Kung)
  * added the fmprb_mat module for matrices of balls of real numbers
  * matrix multiplication
  * interval-aware LU decomposition, solving, inverse and determinant
  * many helper functions and small bugfixes

* 2012-09-14 - version 0.1
* 2012-08-05 - began simplified rewrite
* 2012-04-05 - experimental ball and polynomial code

