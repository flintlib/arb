.. _overview:

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

Arb version 1.1 contains:

* A module (:ref:`fmpr <fmpr>`) for correctly rounded arbitrary-precision
  floating-point arithmetic. Arb numbers have a few special features, such
  as arbitrary-size exponents (useful for combinatorics and asymptotics) and
  dynamic allocation (facilitating implementation of hybrid
  integer/floating-point and mixed-precision algorithms).

* A module (:ref:`fmprb <fmprb>`) for real ball arithmetic, where a ball is
  implemented as a pair of fmpr numbers, and a corresponding module
  (:ref:`fmpcb <fmpcb>`) for complex numbers in rectangular form.

* Functions for fast high-precision evaluation of various
  mathematical constants and special functions, implemented using
  ball arithmetic with rigorous error bounds.

* Modules (:ref:`fmprb_poly <fmprb-poly>`, :ref:`fmpcb_poly <fmpcb-poly>`)
  for polynomials or power series over the real and complex numbers,
  implemented using balls as coefficients,
  with asymptotically fast polynomial multiplication and
  many other operations.

* Modules (:ref:`fmprb_mat <fmprb-mat>`, :ref:`fmpcb_mat <fmpcb-mat>`)
  for matrices over the real and complex numbers,
  implemented using balls as coefficients.
  At the moment, only rudimentary linear algebra operations are provided.

Arb 2.0 adds a new set of types designed for higher performance.
They can be used as drop-in replacements for the types they replace
(minor adjustments may be necessary).

* :ref:`arf <arf>` - replaces :ref:`fmpr <fmpr>` for fixed-precision radii

* :ref:`mag <mag>` - replaces :ref:`fmpr <fmpr>` for arbitrary-precision midpoints

* :ref:`arb <arb>` - replaces :ref:`fmprb <fmprb>` for real numbers

* :ref:`arb_poly <arb-poly>` - replaces :ref:`fmprb_poly <fmprb-poly>` for real polynomials

* :ref:`arb_mat <arb-mat>` - replaces :ref:`fmprb_mat <fmprb-mat>` for real matrices

* :ref:`acb <acb>` - replaces :ref:`fmpcb <fmpcb>` for complex numbers

* :ref:`acb_poly <acb-poly>` - replaces :ref:`fmpcb_poly <fmpcb-poly>` for complex polynomials

* :ref:`acb_mat <acb-mat>` - replaces :ref:`fmpcb <fmpcb>` for complex matrices

Planned features include more transcendental functions and more extensive
polynomial and matrix functionality, as well as further optimizations.

Arb uses `GMP <http://mpir.org>`_ / `MPIR <http://mpir.org>`_ and
`FLINT <http://flintlib.org/>`_
for the underlying integer arithmetic and other functions.
The code conventions borrow from FLINT, and the project might get
merged back into FLINT when the code stabilizes in the future.
Arb also uses `MPFR <http://mpfr.org/>`_ for testing purposes
and for evaluation of some functions.

The current version of Arb implements most of its floating-point arithmetic
naively using high-level FLINT types. The speed at low precision is far from
optimal, and the memory management can sometimes be wasteful. The internals
will be rewritten in the future to fix the inefficiencies,
which eventually should make Arb ball arithmetic about as fast as
mpz or mpfr arithmetic at any precision.

