.. _setup:

Setup
===============================================================================

Download
-------------------------------------------------------------------------------

Tarballs of released versions can be downloaded from https://github.com/fredrik-johansson/arb/releases

Alternatively, you can simply install Arb from a git checkout of https://github.com/fredrik-johansson/arb/.
The master branch is generally safe to use (the test suite should pass at all
times), and recommended for
keeping up with the latest changes.

Dependencies
-------------------------------------------------------------------------------

Arb has the following dependencies:

* Either MPIR (http://www.mpir.org) 2.6.0 or later, or GMP (http://www.gmplib.org) 5.1.0 or later.
  If MPIR is used instead of GMP, it must be compiled with the ``--enable-gmpcompat`` option.
* MPFR (http://www.mpfr.org) 3.0.0 or later.
* FLINT (http://www.flintlib.org) version 2.4 or later. You may also
  use a git checkout of https://github.com/fredrik-johansson/flint2


Installation as part of FLINT
-------------------------------------------------------------------------------

With a sufficiently new version of FLINT, Arb can be compiled as a FLINT
extension package.

Simply put the Arb source directory somewhere, say ``/path/to/arb``.
Then go to the FLINT source directory and build FLINT using::

    ./configure --extensions=/path/to/arb <other options>
    make
    make check       (optional)
    make install

This is convenient, as Arb does not need to be
configured or linked separately. Arb becomes part of the compiled FLINT
library, and the Arb header files will be installed along with the other
FLINT header files.

Standalone installation
-------------------------------------------------------------------------------

To compile, test and install Arb from source as a standalone library,
first install FLINT. Then go to the Arb source directory and run::

    ./configure <options>
    make
    make check       (optional)
    make install

If GMP/MPIR, MPFR or FLINT is installed in some other location than
the default path ``/usr/local``, pass
``--with-gmp=...``, ``--with-mpfr=...`` or ``--with-flint=...`` with
the correct path to configure (type ``./configure --help`` to show
more options).


Running code
-------------------------------------------------------------------------------

Here is an example program to get started using Arb:

.. code-block:: c

    #include "arb.h"

    int main()
    {
        arb_t x;
        arb_init(x);
        arb_const_pi(x, 50 * 3.33);
        arb_printn(x, 50); printf("\n");
        printf("Computed with arb-%s\n", arb_version);
        arb_clear(x);
    }

Compile it with::

    gcc -larb test.c

or (if Arb is built as part of FLINT)::

    gcc -lflint test.c

If the Arb/FLINT header and library files are not in a standard location
(``/usr/local`` on most systems), you may also have to pass options such as::

    -I/path/to/arb -I/path/to/flint -L/path/to/flint -L/path/to/arb

to ``gcc``. Finally, to run the program, make sure that the linker
can find the FLINT (and Arb) libraries. If they are installed in a
nonstandard location, you can for example add this path to the
``LD_LIBRARY_PATH`` environment variable.

The output of the example program should be something like the following::

    [3.1415926535897932384626433832795028841971693993751 +/- 6.28e-50]
    Computed with arb-2.4.0

