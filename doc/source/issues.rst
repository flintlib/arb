.. _issues:

Potential issues
===============================================================================

Interface changes
-------------------------------------------------------------------------------

As this is an early version, note that any part of the interface is
subject to change without warning! Most of the core interface should
be stable at this point, but no guarantees are made.

Correctness
-------------------------------------------------------------------------------

Except where otherwise specified, Arb is designed to produce
provably correct error bounds. The code has been written carefully,
and the library is extensively tested.
However, like any complex mathematical software, Arb is virtually certain to
contains bugs, so the usual precautions are advised:

* Perform sanity checks on the output (check known mathematical relations; recompute to another precision and compare)
* Compare against other mathematical software
* Read the source code to verify that it does what it is supposed to do

All bug reports are highly welcome!

Integer overflow
-------------------------------------------------------------------------------

Machine-size integers are used for precisions, sizes of integers in
bits, lengths of polynomials, and similar quantities that relate
to sizes in memory. Very few checks are performed to verify that
such quantities do not overflow.
Precisions and lengths exceeding a small fraction
of *LONG_MAX*, say `2^{24} \sim 10^7` on 32-bit systems,
should be regarded as resulting in undefined behavior.
On 64-bit systems this should generally not be an issue,
since most calculations will exhaust the available memory
(or the user's patience waiting for the computation to complete)
long before running into integer overflows.
However, the user needs to be wary of unintentionally passing input
parameters of order *LONG_MAX* or negative parameters where
positive parameters are expected, for example due to a runaway loop
that repeatedly increases the precision.

This caveat does not apply to exponents of floating-point numbers,
which are represented as arbitrary-precision integers, nor to
integers used as numerical scalars (e.g. :func:`fmprb_mul_si`).
However, it still applies to conversions and operations where
the result is requested exactly and sizes become an issue.
For example, trying to convert
the floating-point number `2^{2^{100}}` to an integer could
result in anything from a silent wrong value to thrashing followed
by a crash, and it is the user's responsibility not
to attempt such a thing.

Thread safety and caches
-------------------------------------------------------------------------------

Arb should be fully threadsafe, provided that both MPFR and FLINT have
been built in threadsafe mode. Please note that thread safety is
not currently tested, and extra caution when developing
multithreaded code is therefore recommended.

Arb may cache some data (such as the value of `\pi` and
Bernoulli numbers) to speed up various computations. In threadsafe mode,
caches use thread-local storage (there is currently no way to save memory
and avoid recomputation by having several threads share the same cache).
Caches can be freed by calling the ``flint_cleanup()`` function. To avoid
memory leaks, the user should call ``flint_cleanup()`` when exiting a thread.
It is also recommended to call ``flint_cleanup()`` when exiting the main
program (this should result in a clean output when running
`Valgrind <http://valgrind.org/>`_, and can help catching memory issues).

