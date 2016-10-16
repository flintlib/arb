Arb - a C library for arbitrary-precision interval arithmetic
=============================================================

.. only:: latex

    Introduction
    ::::::::::::

    Welcome to Arb's documentation!
    Arb is a C library for arbitrary-precision interval arithmetic,
    using a midpoint-radius representation ("ball arithmetic").
    It supports real and complex numbers, polynomials, power series,
    matrices, and evaluation of many transcendental functions.
    All operations are done with automatic, rigorous error bounds.
    The code is thread-safe, portable, and extensively tested.

    Arb is free software distributed under the
    GNU Lesser General Public License (LGPL), version 2.1 or later
    (see :ref:`license`).

    The git repository is https://github.com/fredrik-johansson/arb/

    Arb is developed by `Fredrik Johansson <http://fredrikj.net/>`_
    (fredrik.johansson@gmail.com), with help from many
    contributors (see :ref:`credits`).
    Questions and discussion about Arb are welcome on the
    `flint-devel <https://groups.google.com/d/forum/flint-devel>`_ mailing list.
    There is also an `issue tracker <https://github.com/fredrik-johansson/arb/issues>`_
    for bug reports and feature requests.
    Development progress is sometimes covered on
    `Fredrik's blog <http://fredrikj.net/blog/>`_.

    This documentation is available in HTML format at http://arblib.org and in
    PDF format at http://arblib.org/arb.pdf.
    The version of the documentation you are currently reading was updated
    |today| and describes Arb |version|.
    Documentation for :ref:`specific release versions <history>`
    is also available in PDF format.

.. only:: html

        .. image:: _static/arbtext.png

    Welcome to Arb's documentation!
    Arb is a C library for arbitrary-precision interval arithmetic,
    using a midpoint-radius representation ("ball arithmetic").
    It supports real and complex numbers, polynomials, power series,
    matrices, and evaluation of many transcendental functions.
    All operations are done with automatic, rigorous error bounds.
    The code is thread-safe, portable, and extensively tested.

    Arb is free software distributed under the
    GNU Lesser General Public License (LGPL), version 2.1 or later
    (see :ref:`license`).

    The git repository is https://github.com/fredrik-johansson/arb/

    Arb is developed by `Fredrik Johansson <http://fredrikj.net/>`_
    (fredrik.johansson@gmail.com), with help from many
    contributors (see :ref:`credits`).
    Questions and discussion about Arb are welcome on the
    `flint-devel <https://groups.google.com/d/forum/flint-devel>`_ mailing list.
    There is also an `issue tracker <https://github.com/fredrik-johansson/arb/issues>`_
    for bug reports and feature requests.
    Development progress is sometimes covered on
    `Fredrik's blog <http://fredrikj.net/blog/>`_.

    This documentation is available in HTML format at http://arblib.org and in
    PDF format at http://arblib.org/arb.pdf.
    The version of the documentation you are currently reading was updated
    |today| and describes Arb |version|.
    Documentation for :ref:`specific release versions <history>`
    is also available in PDF format.

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
   hurwitz.rst
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

