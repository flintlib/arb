Feature overview
===============================================================================

Ball arithmetic, also known as mid-rad interval arithmetic, is an
extension of floating-point arithmetic in which an error bound is
attached to each variable. This allows doing rigorous computations
over the real numbers, while avoiding the overhead of
traditional (inf-sup) interval arithmetic at high precision,
and eliminating much of the need for time-consuming
and bug-prone manual error analysis associated with
standard floating-point arithmetic. (See for example [Hoe2009]_.)

At the moment, Arb contains:

* A module (fmpr) for correctly rounded arbitrary-precision
  floating-point arithmetic. Arb numbers have a few special features, such
  as arbitrary-size exponents (useful for combinatorics and asymptotics) and
  dynamic allocation (facilitating implementation of hybrid
  integer/floating-point and mixed-precision algorithms).

* A module (fmprb) for real ball arithmetic, where a ball is
  implemented as a pair of fmpr numbers.

* Functions for fast high-precision computation of some mathematical constants,
  based on ball arithmetic.

* A module (fmprb_poly) for polynomials or power series over the real numbers,
  implemented using balls as coefficients, with fast polynomial multiplication.

* A rudimentary module (fmprb_mat) for matrices over the real numbers,
  implemented using balls as coefficients.

Planned features include: transcendental functions and more extensive
polynomial and matrix functionality, as well as support for complex numbers.

Arb uses `MPIR <http://mpir.org>`_ and `FLINT <http://flintlib.org/>`_
for the underlying integer arithmetic and other functions.
The code conventions borrow from FLINT, and the project might get
merged back into FLINT when the code stabilizes in the future.
Arb also uses `MPFR <http://mpfr.org/>`_, mainly for testing purposes
and fallback code.

The current version of Arb implements most of its floating-point arithmetic
naively using high-level FLINT types. The speed at low precision is far from
optimal, and the memory management can sometimes be wasteful. The internals
will be rewritten in the future to fix the inefficiencies,
which eventually should make Arb ball arithmetic about as fast as mpz or mpfr arithmetic at any precision.

**Warning**: as this is an early version, any part of the interface is
subject to change! Also be aware that there are known and unknown bugs.
