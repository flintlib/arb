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

Other implementations of ball arithmetic include
`iRRAM <http://irram.uni-trier.de/>`_ and
`Mathemagix <http://www.mathemagix.org/www/mmdoc/doc/html/main/index.en.html>`_.
In contrast to those systems, Arb is more focused on low-level arithmetic
and computation of transcendental functions needed for
number theory. Arb also differs in some technical aspects of
the implementation.

Arb 2.x contains:

* A module (:ref:`arf <arf>`) for correctly rounded arbitrary-precision
  floating-point arithmetic. Arb's floating-point numbers have a few special
  features, such as arbitrary-size exponents (useful for combinatorics and
  asymptotics) and dynamic allocation (facilitating implementation of hybrid
  integer/floating-point and mixed-precision algorithms).

* A module (:ref:`mag <mag>`) for representing magnitudes (error bounds)
  more efficiently than with an arbitrary-precision floating-point type.

* A module (:ref:`arb <arb>`) for real ball arithmetic, where a ball is
  implemented as an *arf* midpoint and a *mag* radius.

* A module (:ref:`acb <acb>`) for complex numbers in rectangular form,
  represented as pairs real balls.

* Functions for fast high-precision evaluation of various
  mathematical constants and special functions, implemented using
  ball arithmetic with rigorous error bounds.

* Modules (:ref:`arb_poly <arb-poly>`, :ref:`acb_poly <acb-poly>`)
  for polynomials or power series over the real and complex numbers,
  implemented using balls as coefficients,
  with asymptotically fast polynomial multiplication and
  many other operations.

* Modules (:ref:`arb_mat <arb-mat>`, :ref:`acb_mat <acb-mat>`)
  for matrices over the real and complex numbers,
  implemented using balls as coefficients.
  At the moment, only rudimentary linear algebra operations are provided.

Arb 1.x used a different set of numerical base types (*fmpr*, *fmprb*
and *fmpcb*). These types had a slightly simpler internal representation,
but generally had worse performance. Almost all methods for the Arb 1.x types
have now been ported to faster equivalents for the Arb 2.x types.
The last version to include both the Arb 1.x and Arb 2.x types and methods
was Arb 2.2. As of Arb 2.3, only a small set of *fmpr* and *fmprb*
methods are left for fallback and testing purposes.

Planned features include more transcendental functions and more extensive
polynomial and matrix functionality, as well as further optimizations.

Arb uses `GMP <http://mpir.org>`_ / `MPIR <http://mpir.org>`_ and
`FLINT <http://flintlib.org/>`_
for the underlying integer arithmetic and other functions.
The code conventions borrow from FLINT, and the project might get
merged back into FLINT when the code stabilizes in the future.
Arb also uses `MPFR <http://mpfr.org/>`_ for testing purposes
and for evaluation of some functions.

