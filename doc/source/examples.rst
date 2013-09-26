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
working precision of roughly `n \log_2(10)` bits.

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

keiper_li.c
-------------------------------------------------------------------------------

Given an input integer *n*, this program rigorously computes numerical
values of the Keiper-Li coefficients
`\lambda_0, \ldots, \lambda_n`. The Keiper-Li coefficients
have the property that `\lambda_n > 0` for all `n > 0` if and only if the
Riemann hypothesis is true. This program was used for the record
computations described in [Joh2013]_ (the paper describes
the algorithm in some more detail).

The program takes the following parameters::

    keiper_li n [-prec prec] [-threads num_threads] [-out out_file]

The program prints the first and last few coefficients. It can optionally
write all the computed data to a file. The working precision defaults
to a value that should give all the coefficients to a few digits of
accuracy, but can optionally be set higher (or lower).
On a multicore system, using several threads results in faster
execution.

Sample output::

    fredrik@lemur:~/src/arb$ build/examples/keiper_li 1000 -threads 2
    zeta: cpu/wall(s): 0.4 0.244
    virt/peak/res/peak(MB): 167.98 294.69 5.09 7.43
    log: cpu/wall(s): 0.03 0.038
    gamma: cpu/wall(s): 0.02 0.016
    binomial transform: cpu/wall(s): 0.01 0.018
    0: -0.69314718055994530941723212145817656807550013436026 +/- 6.5389e-347
    1: 0.023095708966121033814310247906495291621932127152051 +/- 2.0924e-345
    2: 0.046172867614023335192864243096033943387066108314123 +/- 1.674e-344
    3: 0.0692129735181082679304973488726010689942120263932 +/- 5.0219e-344
    4: 0.092197619873060409647627872409439018065541673490213 +/- 2.0089e-343
    5: 0.11510854289223549048622128109857276671349132303596 +/- 1.0044e-342
    6: 0.13792766871372988290416713700341666356138966078654 +/- 6.0264e-342
    7: 0.16063715965299421294040287257385366292282442046163 +/- 2.1092e-341
    8: 0.18321945964338257908193931774721859848998098273432 +/- 8.4368e-341
    9: 0.20565733870917046170289387421343304741236553410044 +/- 7.5931e-340
    10: 0.22793393631931577436930340573684453380748385942738 +/- 7.5931e-339
    991: 2.3196617961613367928373899656994682562101430813341 +/- 2.461e-11
    992: 2.3203766239254884035349896518332550233162909717288 +/- 9.5363e-11
    993: 2.321092061239733282811659116333262802034375592414 +/- 1.8495e-10
    994: 2.3218073540188462110258826121503870112747188888893 +/- 3.5907e-10
    995: 2.3225217392815185726928702951225314023773358152533 +/- 6.978e-10
    996: 2.3232344485814623873333223609413703912358283071281 +/- 1.3574e-09
    997: 2.3239447114886014522889542667580382034526509232475 +/- 2.6433e-09
    998: 2.3246517591032700808344143240352605148856869322209 +/- 5.1524e-09
    999: 2.3253548275861382119812576052060526988544993162101 +/- 1.0053e-08
    1000: 2.3260531616864664574065046940832238158044982041872 +/- 3.927e-08
    virt/peak/res/peak(MB): 170.18 294.69 7.51 7.51


real_roots.c
-------------------------------------------------------------------------------

This program isolates the roots of a function on the interval `(a,b)`
(where *a* and *b* are input as double-precision literals)
using the routines in the :ref:`fmprb_calc <fmprb-calc>` module.
The program takes the following arguments::

    real_roots function a b [-refine d] [-verbose] [-maxdepth n] [-maxeval n] [-maxfound n] [-prec n]

The following functions (specified by an integer code) are implemented:

* 0 - `Z(x)` (Riemann-Siegel Z-function)
* 1 - `\sin(x)`
* 2 - `\sin(x^2)`
* 3 - `\sin(1/x)`

The following options are available:

* ``-refine d``: If provided, after isolating the roots, attempt to refine
  the roots to *d* digits of accuracy using a few bisection steps followed
  by Newton's method with adaptive precision, and then print them.

* ``-verbose``: Print more information.

* ``-maxdepth n``: Stop searching after *n* recursive subdivisions.

* ``-maxeval n``: Stop searching after approximately *n* function evaluations
   (the actual number evaluations will be a small multiple of this).

* ``-maxfound n``: Stop searching after having found *n* isolated roots.

* ``-prec n``: Working precision to use for the root isolation.

With *function* 0, the program isolates roots of the Riemann zeta function
on the critical line, and guarantees that no roots are missed
(there are more efficient ways to do this, but it is a nice example)::

    popeye:/scratch/fjohanss/64/arb> build/examples/real_roots 0 0.0 50.0 -verbose
    interval: 25 +/- 25
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    found isolated root in: 14.12353515625 +/- 0.012207
    found isolated root in: 21.0205078125 +/- 0.024414
    found isolated root in: 25.0244140625 +/- 0.024414
    found isolated root in: 30.43212890625 +/- 0.012207
    found isolated root in: 32.9345703125 +/- 0.024414
    found isolated root in: 37.5732421875 +/- 0.024414
    found isolated root in: 40.9423828125 +/- 0.024414
    found isolated root in: 43.32275390625 +/- 0.012207
    found isolated root in: 48.01025390625 +/- 0.012207
    found isolated root in: 49.76806640625 +/- 0.012207
    ---------------------------------------------------------------
    Found roots: 10
    Subintervals possibly containing undetected roots: 0
    Function evaluations: 3425
    cpu/wall(s): 1.22 1.229
    virt/peak/res/peak(MB): 20.63 20.66 2.23 2.23

Find just one root and refine it to approximately 75 digits::

    popeye:/scratch/fjohanss/64/arb> build/examples/real_roots 0 0.0 50.0 -maxfound 1 -refine 75
    interval: 25 +/- 25
    maxdepth = 30, maxeval = 100000, maxfound = 1, low_prec = 30
    refined root:
    14.134725141734693790457251983562470270784257115699243175685567460149963429809 +/- 8.4532e-81

    ---------------------------------------------------------------
    Found roots: 1
    Subintervals possibly containing undetected roots: 8
    Function evaluations: 992
    cpu/wall(s): 0.41 0.415
    virt/peak/res/peak(MB): 20.76 20.76 2.23 2.23

Find roots of `\sin(x^2)` on `(0,50)`. The algorithm cannot isolate
the root at `x = 0` (it is at the endpoint of the interval, and in any
case a root of multiplicity higher than one). The failure is reported::

    popeye:/scratch/fjohanss/64/arb> build/examples/real_roots 2 0 100
    interval: 50 +/- 50
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 3183
    Subintervals possibly containing undetected roots: 1
    Function evaluations: 34058
    cpu/wall(s): 0.26 0.263
    virt/peak/res/peak(MB): 20.73 20.76 1.72 1.72

This does not miss any roots::

    popeye:/scratch/fjohanss/64/arb> build/examples/real_roots 2 1 100
    interval: 50.5 +/- 49.5
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 3183
    Subintervals possibly containing undetected roots: 0
    Function evaluations: 34039
    cpu/wall(s): 0.26 0.266
    virt/peak/res/peak(MB): 20.73 20.76 1.70 1.70

Looking for roots of `\sin(1/x)` on `(0,1)`, the algorithm finds many roots,
but will never find all of them since there are infinitely many::

    popeye:/scratch/fjohanss/64/arb> build/examples/real_roots 3 0.0 1.0
    interval: 0.5 +/- 0.5
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 10198
    Subintervals possibly containing undetected roots: 24695
    Function evaluations: 202587
    cpu/wall(s): 1.73 1.731
    virt/peak/res/peak(MB): 21.84 22.89 2.76 2.76

Remark: the program always computes rigorous containing intervals
for the roots, but the accuracy after refinement could be less than *d* digits.

