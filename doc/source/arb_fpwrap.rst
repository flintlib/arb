.. _arb_fpwrap:

**arb_fpwrap.h** -- floating-point wrappers of Arb mathematical functions
=========================================================================================

This module provides wrappers of Arb functions intended users who
want accurate floating-point mathematical functions
without necessarily caring about ball arithmetic.
The wrappers take floating-point input, give floating-point output,
and automatically increase the internal working precision
to ensure that the output is accurate
(in the rare case of failure, they output NaN along with an error code).

Supported types:

* ``double`` and ``complex_double`` (53-bit precision)

Limitations:

* The wrappers currently only handle finite input and points where function
  value is finite. For example,
  they do not know that `\log(0) = -\infty` or that `\exp(-\infty) = 0`.
  Singular input or output result in ``FPWRAP_UNABLE`` and a NaN output value.
  Evaluation of limit values may be implemented in the future for some functions.
* The wrappers currently treat ``-0.0`` as ``+0.0``. Users who need to
  distinguish signs of zero, e.g. on branch cuts, currently need to do so
  manually.
* When requesting *correct rounding*, the wrappers can fail to converge
  in asymptotic or exact cases (where special algorithms are required).
* If the value is computed accurately internally but is too small to represent
  as a floating-point number, the result will be ``-0.0`` or ``+0.0`` (on underflow)
  or ``-Inf`` or ``+Inf`` (on overflow). Since the underflowed or overflowed
  result is the best possible floating-point approximation of the true value,
  this outcome is considered correct and the flag ``FPWRAP_SUCCESS`` is returned.
  In the future, return status flags may be added to indicate that underflow
  or overflow has occurred.
* Different rounding modes are not yet implemented.

Option and return flags
-------------------------------------------------------------------------------

Functions return an ``int`` flag indicating the status.

.. macro:: FPWRAP_SUCCESS

    Indicates an accurate result. (Up to inevitable underflow or
    overflow in the final conversion to a floating-point result; see above.)

    This flag has the numerical value 0.

.. macro:: FPWRAP_UNABLE

    Indicates failure (unable to achieve to target accuracy,
    possibly because of a singularity). The output is set to NaN.

    This flag has the numerical value 1.

Functions take a *flags* parameter specifying optional rounding and termination
behavior. This can be set to 0 to use defaults.

.. macro:: FPWRAP_ACCURATE_PARTS

    For complex output, compute both real and imaginary parts to full relative accuracy.
    By default (if this flag is not set), complex results are computed to
    at least 53-bit accuracy as a whole, but if either the real or imaginary
    part is much smaller than the other, that part can have a large relative error.
    Setting this flag can result in slower evaluation or failure to converge
    in some cases.

    This flag has the numerical value 1.

.. macro:: FPWRAP_CORRECT_ROUNDING

    Guarantees *correct rounding*.
    By default (if this flag is not set), real results are accurate up to the
    rounding of the last bit, but the last bit is not guaranteed to
    be rounded optimally.
    Setting this flag can result in slower
    evaluation or failure to converge in some cases.
    Correct rounding automatically applies to both real and imaginary parts
    of complex numbers, so it is unnecessary to set both this flag and
    *FPWRAP_ACCURATE_PARTS*.

    This flag has the numerical value 2.

.. macro:: FPWRAP_WORK_LIMIT

    Multiplied by an integer, specifies the maximum working precision to use
    before giving up. With ``n * FPWRAP_WORK_LIMIT`` added to *flags*,
    `n` levels of precision will be used.
    The default `n = 0` is equivalent to `n = 8`, which for ``double``
    means trying with a working precision of 64, 128, 256, 512, 1024, 2048,
    4096, 8192 bits.
    With ``flags = 2 * FPWRAP_WORK_LIMIT``, we only try 64 and 128
    bits, and with ``flags = 16 * FPWRAP_WORK_LIMIT`` we
    go up to 2097152 bits.

    This flag has the numerical value 65536.

Types
-------------------------------------------------------------------------------

Outputs are passed by reference so that we can return status
flags and so that the interface is uniform for functions with
multiple outputs.

.. type:: complex_double

    A struct of two ``double`` components (``real`` and ``imag``), used to
    represent a machine-precision complex number. We use this custom type
    instead of the complex types defined in ``<complex.h>`` since Arb
    does not depend on C99. Users should easily be able to convert
    to the C99 complex type since the layout in memory is identical.

Functions
-------------------------------------------------------------------------------

.. function:: int arb_fpwrap_double_gamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_gamma(complex_double * res, complex_double x, int flags)

    Gamma function.

.. function:: int arb_fpwrap_double_rgamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_rgamma(complex_double * res, complex_double x, int flags)

    Reciprocal gamma function.

.. function:: int arb_fpwrap_double_lgamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_lgamma(complex_double * res, complex_double x, int flags)

    Log-gamma function.

.. function:: int arb_fpwrap_double_digamma(double * res, double x, int flags)
              int arb_fpwrap_cdouble_digamma(complex_double * res, complex_double x, int flags)

    Digamma function.

.. function:: int arb_fpwrap_double_zeta(double * res, double x, int flags)
              int arb_fpwrap_cdouble_zeta(complex_double * res, complex_double x, int flags)

    Riemann zeta function.

Interfacing from Python
-------------------------------------------------------------------------------

This illustrates how to call functions from Python using ``ctypes``::

    import ctypes
    import ctypes.util

    libarb_path = ctypes.util.find_library('arb')
    libarb = ctypes.CDLL(libarb_path)

    class _complex_double(ctypes.Structure):
        _fields_ = [('real', ctypes.c_double),
                    ('imag', ctypes.c_double)]

    def wrap_double_fun(fun):
        def f(x):
            y = ctypes.c_double()
            if fun(ctypes.byref(y), ctypes.c_double(x), 0):
                raise ValueError(f"unable to evaluate function accurately at {x}")
            return y.value
        return f

    def wrap_cdouble_fun(fun):
        def f(x):
            x = complex(x)
            cx = _complex_double()
            cy = _complex_double()
            cx.real = x.real
            cx.imag = x.imag
            if fun(ctypes.byref(cy), cx, 0):
                raise ValueError(f"unable to evaluate function accurately at {x}")
            return complex(cy.real, cy.imag)
        return f

    zeta = wrap_double_fun(libarb.arb_fpwrap_double_zeta)
    czeta = wrap_cdouble_fun(libarb.arb_fpwrap_cdouble_zeta)

    print(zeta(2.0))
    print(czeta(0.5+1e9j))
    print(zeta(1.0))       # pole, where wrapper throws exception

This should print::

    1.6449340668482264
    (-2.761748029838061-1.6775122409894598j)
    Traceback (most recent call last):
      ...
    ValueError: unable to evaluate function accurately at 1.0

Interfacing from Julia
-------------------------------------------------------------------------------

This illustrates how to call functions from Julia using ``ccall``::

    using Libdl

    dlopen("/home/fredrik/src/arb/libarb.so")

    function zeta(x::Float64)
        cy = Ref{Float64}()
        if Bool(ccall((:arb_fpwrap_double_zeta, :libarb), Cint, (Ptr{Float64}, Float64), cy, x))
            error("unable to evaluate accurately at ", x)
        end
        return cy[]
    end

    function zeta(x::Complex{Float64})
        cy = Ref{Complex{Float64}}()
        if Bool(ccall((:arb_fpwrap_cdouble_zeta, :libarb), Cint, (Ptr{Complex{Float64}}, Complex{Float64}), cy, x))
            error("unable to evaluate accurately at ", x)
        end
        return cy[]
    end

    println(zeta(2.0))
    println(zeta(0.5 + 1e9im))
    println(zeta(1.0))       # pole, where wrapper throws exception

This should print::

    1.6449340668482264
    -2.761748029838061 - 1.6775122409894598im
    ERROR: unable to evaluate accurately at 1.0
    Stacktrace:
     ...

