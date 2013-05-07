Installation and usage basics
===============================================================================

Arb has the following dependencies:

* MPIR (http://www.mpir.org) or GMP (http://www.gmplib.org)
* MPFR (http://www.mpfr.org)
* FLINT (http://www.flintlib.org)

If MPIR is used instead of GMP, it must be compiled with
the --enable-gmpcompat option.

Currently a source checkout of FLINT from
https://github.com/fredrik-johansson/flint2 is required.

To compile, test and install Arb from source, do::

    ./configure <options>
    make
    make check
    make install

If GMP/MPIR, MPFR or FLINT is installed in some other location than
the default path /usr/local, pass the
flag --with-gmp=... --with-mpfr=... or --with-flint=... with
the correct path to configure (type ./configure --help to show
more options).

Here is a simple sample program to get started using Arb:

.. code-block:: c

    #include "fmprb.h"

    int main()
    {
        fmprb_t x;
        fmprb_init(x);
        fmprb_const_pi(x, 50 * 3.33);
        fmprb_printd(x, 50); printf("\n");
        fmprb_clear(x);
    }

The output should be something like the following::

    3.1415926535897932384626433832795028841971693993751 +/- 4.2764e-50

