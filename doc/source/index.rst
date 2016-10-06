Welcome to Arb's documentation!
===============================

.. only:: latex

    Introduction
    ::::::::::::

    Arb is a C library for arbitrary-precision floating-point ball arithmetic,
    developed by `Fredrik Johansson <http://fredrikj.net/>`_
    (fredrik.johansson@gmail.com).
    It supports real and complex numbers, polynomials, power series, matrices,
    and evaluation of many transcendental functions.
    All is done with automatic, rigorous error bounds.

    Arb is free software distributed under the
    GNU Lesser General Public License (LGPL), version 2.1 or later
    (see :ref:`license` for details).

    The git repository is https://github.com/fredrik-johansson/arb/

    The documentation website is http://fredrikj.net/arb/

.. only:: html

        .. image:: _static/arbtext.png

    Arb is a C library for arbitrary-precision floating-point ball arithmetic,
    developed by `Fredrik Johansson <http://fredrikj.net/>`_
    (fredrik.johansson@gmail.com).
    It supports real and complex numbers, polynomials, power series, matrices,
    and evaluation of many transcendental functions.
    All is done with automatic, rigorous error bounds.

    Arb is free software distributed under the
    GNU Lesser General Public License (LGPL), version 2.1 or later
    (see :ref:`license` for details).

    The git repository is https://github.com/fredrik-johansson/arb/

    The documentation website is http://fredrikj.net/arb/

    A `PDF version <http://fredrikj.net/arb/arb.pdf>`_ of this documentation
    is available (:ref:`older versions <history>`).

General information
::::::::::::::::::::

.. toctree::
   :maxdepth: 2

   overview.rst
   setup.rst
   using.rst
   issues.rst
   examples.rst

Floating-point numbers
::::::::::::::::::::::::::::::::::::

The radius and midpoint of a ball are represented using two specialized
floating-point types.

.. toctree::
   :maxdepth: 2

   mag.rst
   arf.rst

Real and complex numbers
::::::::::::::::::::::::::::::::::::

Real numbers (*arb_t*) are represented as midpoint-radius intervals,
also known as balls. Complex numbers (*acb_t*) are represented in rectangular
form, with balls for the real and imaginary parts.

.. toctree::
   :maxdepth: 2

   arb.rst
   acb.rst

Polynomials and power series
::::::::::::::::::::::::::::::::::::

These modules implement dense univariate polynomials with real and complex
coefficients. Truncated power series are supported via methods acting
on polynomials, without introducing a separate power series type.

.. toctree::
   :maxdepth: 2

   arb_poly.rst
   acb_poly.rst

Matrices
::::::::::::::::::::::::::::::::::::

These modules implement dense matrices with real and complex coefficients.
Rudimentary linear algebra is supported.

.. toctree::
   :maxdepth: 2

   arb_mat.rst
   acb_mat.rst

Higher mathematical functions
::::::::::::::::::::::::::::::::::::

These modules implement mathematical functions with complexity
that goes beyond the basics covered directly in the *arb* and *acb*
modules.

.. toctree::
   :maxdepth: 2

   acb_hypgeom.rst
   arb_hypgeom.rst
   acb_modular.rst
   dirichlet.rst
   acb_dirichlet.rst
   bernoulli.rst
   hypgeom.rst
   partitions.rst

Calculus
::::::::::::::::::::::::::::::::::::

Using ball arithmetic, it is possible to do rigorous root-finding and
integration (among other operations)
with generic functions. This code should be considered experimental.

.. toctree::
   :maxdepth: 2

   arb_calc.rst
   acb_calc.rst

Extra utility modules
::::::::::::::::::::::::::::::::::::

Mainly for internal use.

.. toctree::
   :maxdepth: 1

   fmpz_extras.rst
   bool_mat.rst
   dlog.rst
   fmpr.rst

Supplementary algorithm notes
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Here, we give extra proofs, error bounds, and formulas that would be too
lengthy to reproduce in the documentation for each module.

.. toctree::
   :maxdepth: 1

   formulas.rst
   constants.rst
   gamma.rst
   polylogarithms.rst
   hypergeometric.rst
   agm.rst

History, credits and references
:::::::::::::::::::::::::::::::::

.. toctree::
   :maxdepth: 2

   credits.rst

.. toctree::
   :maxdepth: 1

   history.rst

