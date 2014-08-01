.. _history:

History and changes
===============================================================================

For more details, view the commit log
in the git repository https://github.com/fredrik-johansson/arb

* 2014-08-01 - version 2.2.0

  * added functions for computing polylogarithms and order expansions
    of polylogarithms, with support for real and complex s, z
  * added a missing cast affecting C++ compatibility
  * generalized powsum functions to allow a geometric factor
  * improved powsum functions slightly when the exponent is an integer
  * faster arb_log_ui_from_prev
  * added mag_sqrt and mag_rsqrt functions
  * fixed various minor bugs and added missing tests and documentation entries

* 2014-06-20 - version 2.1.0

  * ported most of the remaining functions to the new arb/acb types,
    including:

    * elementary functions (log, atan, etc.)
    * hypergeometric series summation
    * the gamma function
    * the Riemann zeta function and related functions
    * Bernoulli numbers
    * the partition function
    * the calculus modules (rigorous real root isolation, rigorous numerical integration of complex-valued functions)
    * example programs

  * added several missing utility functions to the arf and mag modules

* 2014-05-27 - version 2.0.0

  * new modules mag, arf, arb, arb_poly, arb_mat, acb, acb_poly,
    acb_mat for higher-performance ball arithmetic

  * poly_roots2 and hilbert_matrix2 example programs

  * vector dot product and norm functions (contributed by Abhinav Baid)

* 2014-05-03 - version 1.1.0

  * faster and more accurate error bounds for polynomial multiplication
    (error bounds are now always as good as with classical multiplication,
    and multiplying high-degree polynomials with approximately equal
    coefficients now has proper quasilinear complexity)

  * faster and much less memory-hungry exponentials at very high precision

  * improved the partition function to support n bigger than a single word,
    and enabled the possibility to use two threads for the computation

  * fixed a bug in floating-point arithmetic that caused a too small bound
    for the rounding error to be reported when the result of an inexact
    operation was rounded up to a power of two (this bug did
    not affect the correctness of ball arithmetic, because operations on
    ball midpoints always round down)

  * minor optimizations to floating-point arithmetic

  * improved argument reduction of the digamma function and short series
    expansions of the rising factorial

  * removed the holonomic module for now, as it did not really do anything
    very useful

* 2013-12-21 - version 1.0.0

  * new example programs directory

    * poly_roots example program
    * real_roots example program
    * pi_digits example program
    * hilbert_matrix example program
    * keiper_li example program

  * new fmprb_calc module for calculus with real functions

    * bisection-based root isolation
    * asymptotically fast Newton root refinement

  * new fmpcb_calc module for calculus with complex functions

    * numerical integration using Taylor series

  * scalar functions

    * simplified fmprb_const_euler using published error bound
    * added fmprb_inv
    * fmprb_trim, fmpcb_trim
    * added fmpcb_rsqrt (complex reciprocal square root)
    * fixed bug in fmprb_sqrtpos with nonfinite input
    * slightly improved fmprb powering code
    * added various functions for bounding fmprs by powers of two
    * added fmpr_is_int

  * polynomials and power series

    * implemented scaling to speed up blockwise multiplication
    * slightly faster basecase power series exponentials
    * improved sin/cos/tan/exp for short power series
    * added complex sqrt_series, rsqrt_series
    * implemented the Riemann-Siegel Z and theta functions for real power series
    * added fmprb_poly_pow_series, fmprb_poly_pow_ui and related methods
    * fmprb/fmpcb_poly_contains_fmpz_poly
    * faster composition by monomials
    * implemented Borel transform and binomial transform for real power series

  * matrices

    * implemented matrix exponentials
    * multithreaded fmprb_mat_mul
    * added matrix infinity norm functions
    * added some more matrix-scalar functions
    * added matrix contains and overlaps methods

  * zeta function evaluation

    * multithreaded power sum evaluation
    * faster parameter selection when computing many derivatives
    * implemented binary splitting to speed up computing many derivatives

  * miscellaneous

    * corrections for C++ compatibility (contributed by Jonathan Bober)
    * several minor bugfixes and test code enhancements

* 2013-08-07 - version 0.7

  * floating-point and ball functions

    * documented, added test code, and fixed bugs for various operations involving a ball containing an infinity or NaN
    * added reciprocal square root functions (fmpr_rsqrt, fmprb_rsqrt) based on mpfr_rec_sqrt
    * faster high-precision division by not computing an explicit remainder
    * slightly faster computation of pi by using new reciprocal square root and division code
    * added an fmpr function for approximate division to speed up certain radius operations
    * added fmpr_set_d for conversion from double
    * allow use of doubles to optionally compute the partition function faster but without an error bound
    * bypass mpfr overflow when computing the exponential function to extremely high precision (approximately 1 billion digits)
    * made fmprb_exp faster for large numbers at extremely high precision by skipping the log(2) removal
    * made fmpcb_lgamma faster at high precision by speeding up the argument reduction branch computation
    * added fmprb_asin, fmprb_acos
    * added various other utility functions to the fmprb module
    * added a function for computing the Glaisher constant
    * optimized evaluation of the Riemann zeta function at high precision

  * polynomials and power series

    * made squaring of polynomials faster than generic multiplication
    * implemented power series reversion (various algorithms) for the fmprb_poly type
    * added many fmprb_poly utility functions (shifting, truncating, setting/getting coefficients, etc.)
    * improved power series division when either operand is short
    * improved power series logarithm when the input is short
    * improved power series exponential to use the basecase algorithm for short input regardless of the output size
    * added power series square root and reciprocal square root
    * added atan, tan, sin, cos, sin_cos, asin, acos fmprb_poly power series functions
    * added Newton iteration macros to simplify various functions
    * added gamma functions of real and complex power series ([fmprb/fmpcb]_poly_[gamma/rgamma/lgamma]_series)
    * added wrappers for computing the Hurwitz zeta function of a power series ([fmprb/fmpcb]_poly_zeta_series)
    * implemented sieving and other optimizations to improve performance for evaluating the zeta function of a short power series
    * improved power series composition when the inner series is linear
    * added many fmpcb_poly versions of nearly all fmprb_poly functions
    * improved speed and stability of series composition/reversion by balancing the power table exponents

  * other

    * added support for freeing all cached data by calling flint_cleanup()
    * introduced fmprb_ptr, fmprb_srcptr, fmpcb_ptr, fmpcb_srcptr typedefs for cleaner function signatures
    * various bug fixes and general cleanup

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

