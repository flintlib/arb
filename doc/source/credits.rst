.. _credits:

Credits and references
===============================================================================

Arb is licensed GNU General Public License version 2, or any later version.

Fredrik Johansson is the main author. The project was started in 2012
as a numerical extension of FLINT, and the initial design was heavily based
on FLINT 2.0 (with particular credit to Bill Hart and Sebastian Pancratz).

From 2012 to July 2014, Fredrik's work on Arb was supported by
Austrian Science Fund FWF Grant Y464-N18 (Fast Computer Algebra
for Special Functions).
During that period, he was a PhD student (and briefly a postdoc) at
RISC, Johannes Kepler University, Linz, supervised by Manuel Kauers.

From September 2014 to October 2015, Fredrik was a postdoc at
INRIA Bordeaux and Institut de Mathématiques de Bordeaux,
in the LFANT project-team headed by Andreas Enge. During that period,
Fredrik's work on Arb was supported
by ERC Starting Grant ANTICS 278537 (Algorithmic Number Theory in
Computer Science) http://cordis.europa.eu/project/rcn/101288_en.html
Since October 2015, Fredrik is a CR2 researcher in the LFANT team,
funded by INRIA.

Contributors
-------------------------------------------------------------------------------

Several people have contributed patches, bug reports, or substantial feedback.
This list is probably incomplete.

* Bill Hart - build system, Windows 64 support, design of FLINT
* Sebastian Pancratz - divide-and-conquer polynomial composition algorithm (taken from FLINT)
* The MPFR development team - Arb includes two-limb multiplication code taken from MPFR
* Jonathan Bober - Dirichlet characters (the code in Arb is based on his Cython implementation), C++ compatibility fixes
* Yuri Matiyasevich - feedback about the zeta function and root-finding code
* Abhinav Baid - dot product and norm functions
* Ondřej Čertík - bug reports, feedback
* Andrew Booker - bug reports, feedback
* Francesco Biscani - C++ compatibility fixes, feedback
* Clemens Heuberger - work on Arb interface in Sage, feedback
* Marc Mezzarobba - work on Arb interface in Sage, bug reports, feedback
* Pascal Molin - feedback
* Ricky Farr - convenience functions, feedback
* Marcello Seri - fix for static builds on OS X
* Tommy Hofmann - matrix transpose, comparison, other utility methods, Julia interface
* Alexander Kobel - documentation and code cleanup patches
* Hrvoje Abraham - patches for MinGW compatibility
* Julien Puydt - soname versioning support
* Alex Griffing - sinc function, matrix trace, improved matrix squaring, boolean matrices, improved structured matrix exponentials, miscellaneous patches
* Jeroen Demeyer - patch for major bug on PPC64

Software
-------------------------------------------------------------------------------

The following software has been helpful in the development of Arb.

* GMP (Torbjörn Granlund and others), http://gmplib.org
* MPIR (Brian Gladman, Jason Moxham, William Hart and others), http://mpir.org
* MPFR (Guillaume Hanrot, Vincent Lefèvre, Patrick Pélissier, Philippe Théveny, Paul Zimmermann and others), http://mpfr.org
* FLINT (William Hart, Sebastian Pancratz, Andy Novocin, Fredrik Johansson, David Harvey and others), http://flintlib.org
* Sage (William Stein and others), http://sagemath.org
* Pari/GP (The Pari group), http://pari.math.u-bordeaux.fr/
* SymPy (Ondřej Čertík, Aaron Meurer and others), http://sympy.org
* mpmath (Fredrik Johansson and others), http://mpmath.org
* Mathematica (Wolfram Research), http://www.wolfram.com/mathematica
* HolonomicFunctions (Christoph Koutschan), http://www.risc.jku.at/research/combinat/software/HolonomicFunctions/
* Sphinx (George Brandl and others), http://sphinx.pocoo.org
* CM (Andreas Enge), http://www.multiprecision.org/index.php?prog=cm

Citing Arb
-------------------------------------------------------------------------------

If you wish to cite Arb in a scientific paper, the following reference can be used (you may also cite the manual or the website directly):

\F. Johansson. "Arb: a C library for ball arithmetic", *ACM Communications in Computer Algebra*, 47(4):166-169, 2013.

In BibTeX format::

  @article{Johansson2013arb,
    title={{A}rb: a {C} library for ball arithmetic},
    author={F. Johansson},
    journal={ACM Communications in Computer Algebra},
    volume={47},
    number={4},
    pages={166--169},
    year={2013},
    publisher={ACM}
  }


Bibliography
-------------------------------------------------------------------------------

(In the PDF edition, this section is empty. See the bibliography listing at the end of the document.)

.. [Arn2010] \J. Arndt, *Matters Computational*, Springer (2010), http://www.jjj.de/fxt/#fxtbook

.. [BBC1997] \D. H. Bailey, J. M. Borwein and R. E. Crandall, "On the Khintchine constant", Mathematics of Computation 66 (1997) 417-431

.. [Blo2009] \R. Bloemen, "Even faster zeta(2n) calculation!", https://web.archive.org/web/20141101133659/http://xn--2-umb.com/09/11/even-faster-zeta-calculation

.. [BBC2000] \J. Borwein, D. M. Bradley and R. E. Crandall, "Computational strategies for the Riemann zeta function", Journal of Computational and Applied Mathematics 121 (2000) 247-296

.. [BZ1992]_ \J. Borwein and I. Zucker, "Fast evaluation of the gamma function for small rational fractions using complete elliptic integrals of the first kind", IMA Journal of Numerical Analysis 12 (1992) 519-526

.. [Bor1987]_ \P. Borwein, "Reduced complexity evaluation of hypergeometric functions", Journal of Approximation Theory 50:3 (1987)

.. [Bor2000] \P. Borwein, "An Efficient Algorithm for the Riemann Zeta Function", Constructive experimental and nonlinear analysis, CMS Conference Proc. 27 (2000) 29-34, http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf

.. [BM1980] \R. P. Brent and E. M. McMillan, "Some new algorithms for high-precision computation of Euler's constant", Mathematics of Computation 34 (1980) 305-312.

.. [Bre1978] \R. P. Brent, "A Fortran multiple-precision arithmetic package", ACM Transactions on Mathematical Software, 4(1):57–70, March 1978.

.. [Bre2010] \R. P. Brent, "Ramanujan and Euler's Constant", http://wwwmaths.anu.edu.au/~brent/pd/Euler_CARMA_10.pdf

.. [BJ2013] \R. P. Brent and F. Johansson, "A bound for the error term in the Brent-McMillan algorithm", preprint (2013), http://arxiv.org/abs/1312.0039

.. [BZ2011] \R. P. Brent and P. Zimmermann, *Modern Computer Arithmetic*, Cambridge University Press (2011), http://www.loria.fr/~zimmerma/mca/pub226.html

.. [CP2005] \R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, second edition, Springer (2005).

.. [Dup2006] \R. Dupont. "Moyenne arithmético-géométrique, suites de Borchardt et applications." These de doctorat, École polytechnique, Palaiseau (2006). http://http://www.lix.polytechnique.fr/Labo/Regis.Dupont/these_soutenance.pdf

.. [EM2004] \O. Espinosa and V. Moll, "A generalized polygamma function", Integral Transforms and Special Functions (2004), 101-115.

.. [Fil1992] \S. Fillebrown, "Faster Computation of Bernoulli Numbers", Journal of Algorithms 13 (1992) 431-445

.. [GG2003] \J. von zur Gathen and J. Gerhard, *Modern Computer Algebra*, second edition, Cambridge University Press (2003)

.. [GS2003] \X. Gourdon and P. Sebah, "Numerical evaluation of the Riemann Zeta-function" (2003), http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf

.. [HZ2004] \G. Hanrot and P. Zimmermann, "Newton Iteration Revisited" (2004), http://www.loria.fr/~zimmerma/papers/fastnewton.ps.gz

.. [Hoe2009] \J. van der Hoeven, "Ball arithmetic", Technical Report, HAL 00432152 (2009), http://www.texmacs.org/joris/ball/ball-abs.html

.. [Joh2012] \F. Johansson, "Efficient implementation of the Hardy-Ramanujan-Rademacher formula", LMS Journal of Computation and Mathematics, Volume 15 (2012), 341-359, http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=8710297

.. [Joh2013] \F. Johansson, "Rigorous high-precision computation of the Hurwitz zeta function and its derivatives", Numerical Algorithms, http://arxiv.org/abs/1309.2877 http://dx.doi.org/10.1007/s11075-014-9893-1

.. [Joh2014a] \F. Johansson, *Fast and rigorous computation of special functions to high precision*, PhD thesis, RISC, Johannes Kepler University, Linz, 2014. http://fredrikj.net/thesis/

.. [Joh2014b] \F. Johansson, "Evaluating parametric holonomic sequences using rectangular splitting", ISSAC 2014, 256-263. http://dx.doi.org/10.1145/2608628.2608629

.. [Joh2014c] \F. Johansson, "Efficient implementation of elementary functions in the medium-precision range", http://arxiv.org/abs/1410.7176

.. [Joh2015] \F. Johansson, "Computing Bell numbers", http://fredrikj.net/blog/2015/08/computing-bell-numbers/

.. [Kar1998] \E. A. Karatsuba, "Fast evaluation of the Hurwitz zeta function and Dirichlet L-series", Problems of Information Transmission 34:4 (1998), 342-353, http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=425&option_lang=eng

.. [Kob2010] \A. Kobel, "Certified Complex Numerical Root Finding", Seminar on Computational Geometry and Geometric Computing (2010), http://www.mpi-inf.mpg.de/departments/d1/teaching/ss10/Seminar_CGGC/Slides/02_Kobel_NRS.pdf

.. [Kri2013] \A. Krishnamoorthy and D. Menon, "Matrix Inversion Using Cholesky Decomposition" Proc. of the International Conference on Signal Processing Algorithms, Architectures, Arrangements, and Applications (SPA-2013), pp. 70-72, 2013.

.. [MPFR2012] The MPFR team, "MPFR Algorithms" (2012), http://www.mpfr.org/algo.html

.. [NIST2012] National Institute of Standards and Technology, *Digital Library of Mathematical Functions* (2012), http://dlmf.nist.gov/

.. [Olv1997] \F. Olver, *Asymptotics and special functions*, AKP Classics, AK Peters Ltd., Wellesley, MA, 1997. Reprint of the 1974 original.

.. [Rad1973] \H. Rademacher, *Topics in analytic number theory*, Springer, 1973.

.. [PS1973] \M. S. Paterson and L. J. Stockmeyer, "On the number of nonscalar multiplications necessary to evaluate polynomials", SIAM J. Comput (1973)

.. [Smi2001] \D. M. Smith, "Algorithm: Fortran 90 Software for Floating-Point Multiple Precision Arithmetic, Gamma and Related Functions", Transactions on Mathematical Software 27 (2001) 377-387, http://myweb.lmu.edu/dmsmith/toms2001.pdf

.. [Tak2000] \D. Takahashi, "A fast algorithm for computing large Fibonacci numbers", Information Processing Letters 75 (2000) 243-246, http://www.ii.uni.wroc.pl/~lorys/IPL/article75-6-1.pdf


