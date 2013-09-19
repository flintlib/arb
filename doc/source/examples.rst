.. _examples:

Example programs
===============================================================================

The *examples* directory
(https://github.com/fredrik-johansson/arb/tree/master/examples)
contains several complete C programs, which are documented below. Running::

    make examples

will compile the programs and place the binaries in ``build/examples``.

pi.c
-------------------------------------------------------------------------------

This program computes `\pi` to an accuracy of roughly *n* decimal digits
by calling the :func:`fmprb_const_pi` function with a
working precision roughly `n \log_2(10)` bits.

Sample output, computing `\pi` to one million digits::

    fredrik@lemur:~/src/arb$ build/examples/pi 1000000
    computing pi with a precision of 3321933 bits... cpu/wall(s): 0.58 0.586
    virt/peak/res/peak(MB): 28.24 36.84 8.86 15.56
    3.141592654 +/- 1.335e-1000001

The program prints a decimal approximation of the computed ball,
with the midpoint rounded to a number of decimal digits that can be
passed as a second parameter to the program (default = 10).
In the present implementation (see :func:`fmprb_printd`), the
digits are not guaranteed to be correctly rounded.

hilbert_matrix.c
-------------------------------------------------------------------------------

Given an input integer *n*, this program accurately computes the
determinant of the *n* by *n* Hilbert matrix.
Hilbert matrices are notoriously ill-conditioned: although the
entries are close to unit magnitude, the determinant `h_n`
decreases superexponentially (nearly as `1/4^{n^2}`) as
a function of *n*.
This program automatically doubles the working precision
until the ball computed for `h_n` by :func:`fmprb_mat_det`
does not contain zero.

Sample output::

    fredrik@lemur:~/src/arb$ build/examples/hilbert_matrix 200
    prec=20: 0 +/- 5.5777e-330
    prec=40: 0 +/- 2.5785e-542
    prec=80: 0 +/- 8.1169e-926
    prec=160: 0 +/- 2.8538e-1924
    prec=320: 0 +/- 6.3868e-4129
    prec=640: 0 +/- 1.7529e-8826
    prec=1280: 0 +/- 1.8545e-17758
    prec=2560: 2.955454297e-23924 +/- 6.4586e-24044
    success!
    cpu/wall(s): 9.06 9.095
    virt/peak/res/peak(MB): 55.52 55.52 35.50 35.50

