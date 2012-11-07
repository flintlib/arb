History and changes
===============================================================================

For more details, view the detailed commit log
in the git repository https://github.com/fredrik-johansson/arb

* version 0.3-git

  * converted documentation to sphinx
  * initial support for complex numbers and polynomials

    * root isolation for complex polynomials

  * new module fmpz_holonomic for functions/sequences
    defined by linear differential/difference equations
    with polynomial coefficients

    * functions for creating various special sequences and functions
    * some closure properties for sequences
    * Taylor series expansion for differential equations
    * computing the nth entry of a sequence using binary splitting
    * computing the nth entry mod p using fast multipoint evaluation

  * generic code for evaluating hypergeometric series
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

