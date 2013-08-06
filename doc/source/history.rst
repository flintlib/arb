.. _history:

History and changes
===============================================================================

For more details, view the detailed commit log
in the git repository https://github.com/fredrik-johansson/arb

* 2013-05-31 - version 0.6

  * made fast polynomial multiplication over the reals numerically stable by using a blockwise algorithm
  * disabled default use of the Gauss formula for multiplication of complex polynomials, to improve numerical stability
  * added division and remainder for complex polynomials
  * added fast multipoint evaluation and interpolation for complex polynomials
  * added missing fmprb_poly_sub and fmpcb_poly_sub functions
  * faster exponentials (fmprb_exp and dependent functions) at low precision, using precomputation
  * rewrote fmpr_add and fmpr_sub using mpn level code, improving efficiency at low precision
  * ported the partition function implementation from flint (using ball arithmetic
    in all steps of the calculation to guarantee correctness)
  * ported algorithm for computing the cosine minimal polynomial from flint (using
    ball arithmetic to guarantee correctness)
  * support using gmp instead of mpir
  * only use thread-local storage when enabled in flint
  * slightly faster error bounding for the zeta function
  * added some other helper functions

* 2013-03-28 - version 0.5

  * arithmetic and elementary functions

    * added fmpr_get_fmpz, fmpr_get_si
    * fixed accuracy problem with fmprb_div_2expm1
    * special-cased squaring of complex numbers
    * added various fmpcb convenience functions (addmul_ui, etc)
    * optimized fmpr_cmp_2exp_si and fmpr_cmpabs_2exp_si, and added test code for comparison functions
    * added fmprb_atan2, also fixing a bug in fmpcb_arg
    * added fmprb_sin_pi, cos_pi, sin_cos_pi etc.
    * added fmprb_sin_pi_fmpq (etc.) using algebraic methods for fast evaluation of roots of unity
    * faster fmprb_poly_evaluate and evaluate_fmpcb using rectangular splitting
    * added fmprb_poly_evaluate2, evaluate2_fmpcb for simultaneously evaluating the derivative
    * added fmprb_poly root polishing code using near-optimal Newton steps (experimental)
    * added fmpr_root, fmprb_root (currently based on MPFR)
    * added fmpr_min, fmpr_max
    * added fmprb_set_interval_fmpr, fmprb_union
    * added fmpr_bits, fmprb_bits, fmpcb_bits for obtaining the mantissa width
    * added fmprb_hypot
    * added complex square roots
    * improved fmprb_log to slightly improve speed, and properly support huge arguments
    * fixed exp, cosh, sinh to work with huge arguments
    * added fmprb_expm1
    * fixed sin, cos, atan to work with huge arguments
    * improved fmprb_pow and fmpcb_pow, including automatic detection of small integer and half-integer exponents
    * added many more elementary functions: fmprb_tan/cot/tanh/coth, fmpcb_tan/cot, and pi versions
    * added fmprb const_e, const_log2, const_log10, const_catalan
    * fixed ball containment/overlap checking to work operate efficiently and correctly with huge exponents
    * strengthened test code for many core operations

  * special functions

    * reorganized zeta function related code
    * faster evaluation of the Riemann zeta function via sieving
    * documented and improved efficiency of the zeta constant binary splitting code
    * calculate error bound in Borwein's algorithm with fmprs instead of using doubles
    * optimized divisions in zeta evaluation via the Euler product
    * use functional equation for Riemann zeta function of a negative argument
    * compute single Bernoulli numbers using ball arithmetic instead of relying on the floating-point code in flint
    * initial code for evaluating the gamma function using its Taylor series
    * much faster rising factorials at high precision, using difference polynomials
    * much faster gamma function at high precision
    * added complex gamma function, log gamma function, and other versions
    * added fmprb_agm (real arithmetic-geometric mean)
    * added fmprb_gamma_fmpq, supporting rapid computation of gamma(p/q) for q = 1,2,3,4,6
    * added real and complex digamma function
    * fixed unnecessary recomputation of Bernoulli numbers
    * optimized computation of Euler's constant, and added proper error bounds
    * avoid reliance on doubles in the hypergeometric series tail bound
    * cleaned up factorials and binomials, computing factorials via gamma

  * other

    * added an fmpz_extras module to collect various internal fmpz helper functions
    * fixed detection of flint header files
    * fixed various other small bugs

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

