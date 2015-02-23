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
by calling the :func:`arb_const_pi` function with a
working precision of roughly `n \log_2(10)` bits.

Sample output, computing `\pi` to one million digits::

    > build/examples/pi 1000000
    computing pi with a precision of 3321933 bits... cpu/wall(s): 0.58 0.586
    virt/peak/res/peak(MB): 28.24 36.84 8.86 15.56
    [3.14159265358979323846{...999959 digits...}42209010610577945815 +/- 3e-1000000]

The program prints an interval guaranteed to contain `\pi`, and where
all displayed digits are correct up to an error of plus or minus
one unit in the last place (see :func:`arb_printn`).
By default, only the first and last few digits are printed.
Pass 0 as a second argument to print all digits (or pass *m* to
print *m* + 1 leading and *m* trailing digits, as above with
the default *m* = 20).

hilbert_matrix.c
-------------------------------------------------------------------------------

Given an input integer *n*, this program accurately computes the
determinant of the *n* by *n* Hilbert matrix.
Hilbert matrices are notoriously ill-conditioned: although the
entries are close to unit magnitude, the determinant `h_n`
decreases superexponentially (nearly as `1/4^{n^2}`) as
a function of *n*.
This program automatically doubles the working precision
until the ball computed for `h_n` by :func:`arb_mat_det`
does not contain zero.

Sample output::

    > build/examples/hilbert_matrix 200
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

    > build/examples/keiper_li 1000 -threads 2
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
using the routines in the :ref:`arb_calc <arb-calc>` module.
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

    > build/examples/real_roots 0 0.0 50.0 -verbose
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

    > build/examples/real_roots 0 0.0 50.0 -maxfound 1 -refine 75
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

Find roots of `\sin(x^2)` on `(0,100)`. The algorithm cannot isolate
the root at `x = 0` (it is at the endpoint of the interval, and in any
case a root of multiplicity higher than one). The failure is reported::

    > build/examples/real_roots 2 0 100
    interval: 50 +/- 50
    maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
    ---------------------------------------------------------------
    Found roots: 3183
    Subintervals possibly containing undetected roots: 1
    Function evaluations: 34058
    cpu/wall(s): 0.26 0.263
    virt/peak/res/peak(MB): 20.73 20.76 1.72 1.72

This does not miss any roots::

    > build/examples/real_roots 2 1 100
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

    > build/examples/real_roots 3 0.0 1.0
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

poly_roots.c
-------------------------------------------------------------------------------

This program finds the complex roots of an integer polynomial
by calling :func:`acb_poly_find_roots` with increasing
precision until the roots certainly have been isolated.
The program takes the following arguments::

    poly_roots [-refine d] [-print d] <poly>

    Isolates all the complex roots of a polynomial with
    integer coefficients. For convergence, the input polynomial
    is required to be squarefree.

    If -refine d is passed, the roots are refined to an absolute
    tolerance better than 10^(-d). By default, the roots are only
    computed to sufficient accuracy to isolate them.
    The refinement is not currently done efficiently.

    If -print d is passed, the computed roots are printed to
    d decimals. By default, the roots are not printed.

    The polynomial can be specified by passing the following as <poly>:

    a <n>          Easy polynomial 1 + 2x + ... + (n+1)x^n
    t <n>          Chebyshev polynomial T_n
    u <n>          Chebyshev polynomial U_n
    p <n>          Legendre polynomial P_n
    c <n>          Cyclotomic polynomial Phi_n
    s <n>          Swinnerton-Dyer polynomial S_n
    b <n>          Bernoulli polynomial B_n
    w <n>          Wilkinson polynomial W_n
    e <n>          Taylor series of exp(x) truncated to degree n
    m <n> <m>      The Mignotte-like polynomial x^n + (100x+1)^m, n > m
    c0 c1 ... cn   c0 + c1 x + ... + cn x^n where all c:s are specified integers

This finds the roots of the Wilkinson polynomial with roots at the
positive integers 1, 2, ..., 100::

    > build/examples/poly_roots -print 15 w 100
    prec=53: 0 isolated roots | cpu/wall(s): 0.42 0.426
    prec=106: 0 isolated roots | cpu/wall(s): 1.37 1.368
    prec=212: 0 isolated roots | cpu/wall(s): 1.48 1.485
    prec=424: 100 isolated roots | cpu/wall(s): 0.61 0.611
    done!
    (1 + 1.7285178043492e-125j)  +/-  (7.2e-122, 7.2e-122j)
    (2 + 5.1605530263601e-122j)  +/-  (3.77e-118, 3.77e-118j)
    (3 + -2.58115555871665e-118j)  +/-  (5.72e-115, 5.72e-115j)
    (4 + 1.02141628524271e-115j)  +/-  (4.38e-112, 4.38e-112j)
    (5 + 1.61326834094948e-113j)  +/-  (2.6e-109, 2.6e-109j)
        ...
    (95 + 4.15294196875447e-62j)  +/-  (6.66e-59, 6.66e-59j)
    (96 + 3.54502401922667e-64j)  +/-  (7.37e-60, 7.37e-60j)
    (97 + -1.67755595325625e-65j)  +/-  (6.4e-61, 6.4e-61j)
    (98 + 2.04638822325299e-65j)  +/-  (4e-62, 4e-62j)
    (99 + -2.73425468028238e-66j)  +/-  (1.71e-63, 1.71e-63j)
    (100 + -1.00950111302288e-68j)  +/-  (3.24e-65, 3.24e-65j)
    cpu/wall(s): 3.88 3.893

This finds the roots of a Bernoulli polynomial which has both real
and complex roots. Note that the program does not attempt to determine
that the imaginary parts of the real roots really are zero (this could
be done by verifying sign changes)::

    > build/examples/poly_roots -refine 100 -print 20 b 16
    prec=53: 16 isolated roots | cpu/wall(s): 0 0.007
    prec=106: 16 isolated roots | cpu/wall(s): 0 0.004
    prec=212: 16 isolated roots | cpu/wall(s): 0 0.004
    prec=424: 16 isolated roots | cpu/wall(s): 0 0.004
    done!
    (-0.94308706466055783383 + -5.512272663168484603e-128j)  +/-  (2.2e-125, 2.2e-125j)
    (-0.75534059252067985752 + 1.937401283040249068e-128j)  +/-  (1.09e-125, 1.09e-125j)
    (-0.24999757119077421009 + -4.5347924422246038692e-130j)  +/-  (3.6e-127, 3.6e-127j)
    (0.24999757152512726002 + 4.2191300761823281708e-129j)  +/-  (4.98e-127, 4.98e-127j)
    (0.75000242847487273998 + 9.0360649917413170142e-128j)  +/-  (8.88e-126, 8.88e-126j)
    (1.2499975711907742101 + 7.8804123808107088267e-127j)  +/-  (2.66e-124, 2.66e-124j)
    (1.7553405925206798575 + 5.432465269253967768e-126j)  +/-  (6.23e-123, 6.23e-123j)
    (1.9430870646605578338 + 3.3035377342500953239e-125j)  +/-  (7.05e-123, 7.05e-123j)
    (-0.99509334829256233279 + 0.44547958157103608805j)  +/-  (5.5e-125, 5.5e-125j)
    (-0.99509334829256233279 + -0.44547958157103608805j)  +/-  (5.46e-125, 5.46e-125j)
    (1.9950933482925623328 + 0.44547958157103608805j)  +/-  (1.44e-122, 1.44e-122j)
    (1.9950933482925623328 + -0.44547958157103608805j)  +/-  (1.43e-122, 1.43e-122j)
    (-0.92177327714429290564 + -1.0954360955079385542j)  +/-  (9.31e-125, 9.31e-125j)
    (-0.92177327714429290564 + 1.0954360955079385542j)  +/-  (1.02e-124, 1.02e-124j)
    (1.9217732771442929056 + 1.0954360955079385542j)  +/-  (9.15e-123, 9.15e-123j)
    (1.9217732771442929056 + -1.0954360955079385542j)  +/-  (8.12e-123, 8.12e-123j)
    cpu/wall(s): 0.02 0.02

