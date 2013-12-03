.. _zeta:

**zeta.h** -- support for the zeta function
===============================================================================

This module implements various algorithms for evaluating the
Riemann zeta function and related functions. The functions provided here are
mainly intended for internal use, though they may be useful to call directly
in some applications where the default algorithm choices are suboptimal.
Most applications should use the user-friendly functions
in the :ref:`fmprb <fmprb>` and :ref:`fmpcb <fmpcb>` modules (or for
power series, the functions in the
:ref:`fmprb_poly <fmprb-poly>` and :ref:`fmpcb_poly <fmpcb-poly>`
modules).

Integer zeta values
-------------------------------------------------------------------------------

.. function:: void zeta_apery_bsplit(fmprb_t x, long prec)

    Sets *x* to Apery's constant `\zeta(3)`, computed by applying binary
    splitting to a hypergeometric series.

.. function:: void zeta_ui_asymp(fmprb_t z, ulong s, long prec)

    Assuming `s \ge 2`, approximates `\zeta(s)` by `1 + 2^{-s}` along with
    a correct error bound. We use the following bounds: for `s > b`,
    `\zeta(s) - 1 < 2^{-b}`, and generally,
    `\zeta(s) - (1 + 2^{-s}) < 2^{2-\lfloor 3 s/2 \rfloor}`.

.. function:: void zeta_ui_euler_product(fmprb_t z, ulong s, long prec)

    Computes `\zeta(s)` using the Euler product. This is fast only if *s*
    is large compared to the precision.

    Writing `P(a,b) = \prod_{a \le p \le b} (1 - p^{-s})`, we have
    `1/\zeta(s) = P(a,M) P(M+1,\infty)`.

    To bound the error caused by truncating
    the product at `M`, we write `P(M+1,\infty) = 1 - \epsilon(s,M)`.
    Since `0 < P(a,M) \le 1`, the absolute error for `\zeta(s)` is
    bounded by `\epsilon(s,M)`.

    According to the analysis in [Fil1992]_, it holds for all `s \ge 6` and `M \ge 1`
    that `1/P(M+1,\infty) - 1 \le f(s,M) \equiv 2 M^{1-s} / (s/2 - 1)`.
    Thus, we have `1/(1-\epsilon(s,M)) - 1 \le f(s,M)`, and expanding
    the geometric series allows us to conclude that
    `\epsilon(M) \le f(s,M)`.

.. function:: void zeta_ui_bernoulli(fmprb_t x, ulong n, long prec)

    Computes `\zeta(n)` for even *n* via the corresponding Bernoulli number.

.. function:: void zeta_ui_vec_borwein(fmprb_ptr z, ulong start, long num, ulong step, long prec)

    Evaluates `\zeta(s)` at `\mathrm{num}` consecutive integers *s* beginning
    with *start* and proceeding in increments of *step*.
    Uses Borwein's formula ([Bor2000]_, [GS2003]_),
    implemented to support fast multi-evaluation
    (but also works well for a single *s*).

    Requires `\mathrm{start} \ge 2`. For efficiency, the largest *s*
    should be at most about as
    large as *prec*. Arguments approaching *LONG_MAX* will cause
    overflows.
    One should therefore only use this function for *s* up to about *prec*, and
    then switch to the Euler product.

    The algorithm for single *s* is basically identical to the one used in MPFR
    (see [MPFR2012]_ for a detailed description).
    In particular, we evaluate the sum backwards to avoid storing more than one
    `d_k` coefficient, and use integer arithmetic throughout since it
    is convenient and the terms turn out to be slightly larger than
    `2^\mathrm{prec}`.
    The only numerical error in the main loop comes from the division by `k^s`,
    which adds less than 1 unit of error per term.
    For fast multi-evaluation, we repeatedly divide by `k^{\mathrm{step}}`.
    Each division reduces the input error and adds at most 1 unit of
    additional rounding error, so by induction, the error per term
    is always smaller than 2 units.

.. function:: void zeta_ui_borwein_bsplit(fmprb_t x, ulong s, long prec)

    Computes `\zeta(s)` for arbitrary `s \ge 2` using a binary splitting
    implementation of Borwein's algorithm. This has quasilinear complexity
    with respect to the precision (assuming that `s` is fixed).
    We have

    .. math ::

        \zeta(s) = \frac{1}{d_n (1-2^{1-s})}
        \sum_{k=0}^{n-1} \frac{(-1)^k(d_n-d_k)}{(k+1)^s} + \gamma_n(s)

    where

    .. math ::

        d_k = n \sum_{i=0}^k \frac{(n+i-1)! 4^i}{(n-i)! (2i)!}.

    On the domain of interest, `|\gamma_n(s)| \le 3 / (3 + \sqrt 8)^n`.

    We write the summation as a system of first-order recurrences for
    `(s_k, d_k, t_k)` where `t_k = d_k - d_{k-1}`. This system is
    described by the matrix equation

    .. math ::

        \begin{pmatrix} s_{k+1} \\ d_{k+2} \\ t_{k+3} \end{pmatrix}
        =
        \begin{pmatrix}
        1 & (-1)^k (k+1)^{-s} & 0 \\
        0 & 1 & 1 \\
        0 & 0 & u(k)
        \end{pmatrix}
        \begin{pmatrix} s_k \\ d_{k+1} \\ t_{k+2} \end{pmatrix}.

    We derive the binary splitting scheme by considering a product
    of an arbitrary pair in the chain `M_{n-1} M_{n-2} \cdots M_1 M_0`.
    This gives

    .. math ::

        \begin{pmatrix}
        1 & A_L & B_L \\
        0 & 1 & C_L \\
        0 & 0 & D_L
        \end{pmatrix}
        \begin{pmatrix}
        1 & A_R & B_R \\
        0 & 1 & C_R \\
        0 & 0 & D_R
        \end{pmatrix} =
        \begin{pmatrix}
        1 & A_L+A_R & B_R+A_L C_R+B_L D_R \\
        0 & 1 & C_R+C_L D_R \\
        0 & 0 & D_L D_R
        \end{pmatrix}.

    The next step is to clear denominators. Instead of putting the
    whole matrix on a common denominator, we optimize by putting `C, D` on a
    denominator `Q_1` (the product of denominators of `u`) and `A, B` on
    a common denominator `Q_3 = Q_1 Q_2` (where `Q_2` is the product of
    `(k+1)^s` factors). This gives a small efficiency improvement. Thus,
    we have

    .. math ::

        \begin{pmatrix}
        1 & \dfrac{A_L}{Q_{3L}} & \dfrac{B_L}{Q_{3L}} \\[3ex]
        0 & 1 & \dfrac{C_L}{Q_{1L}} \\[3ex]
        0 & 0 & \dfrac{D_L}{Q_{1L}}
        \end{pmatrix}
        \begin{pmatrix}
        1 & \dfrac{A_R}{Q_{3R}} & \dfrac{B_R}{Q_{3R}} \\[3ex]
        0 & 1 & \dfrac{C_R}{Q_{1R}} \\[3ex]
        0 & 0 & \dfrac{D_R}{Q_{1R}}
        \end{pmatrix} =
        \begin{pmatrix}
        1 & \dfrac{Q_{3L} A_R + A_L Q_{3R}}{Q_{3L} Q_{3R}} & \dfrac{Q_{3L} B_R + A_L C_R Q_{2R} + B_L D_R Q_{2R}}{Q_{3L} Q_{3R}} \\[3ex]
        0 & 1 & \dfrac{Q_{1L} C_R + C_L D_R}{Q_{1L} Q_{1R}} \\[3ex]
        0 & 0 & \dfrac{D_L D_R}{Q_{1L} Q_{1R}}
        \end{pmatrix}.

    In the final matrix, we note that 
    `A / Q_3 = \sum_k (-1)^k (k+1)^{-s}`, and `C / Q_1 = d_n`.
    Thus `(1 / d_n) \sum_k (-1)^k (k+1)^{-s} (d_n - d_k)` is given by
    `A/Q_3 - (B/Q_3) / (C/Q_1) = (A C - B Q_1) / (Q_3 C)`.

.. function:: void zeta_ui(fmprb_t x, ulong s, long prec)

    Computes `\zeta(s)` for nonnegative integer `s \ne 1`, automatically
    choosing an appropriate algorithm.

.. function:: void zeta_ui_vec(fmprb_ptr x, ulong start, long num, long prec)

.. function:: void zeta_ui_vec_even(fmprb_ptr x, ulong start, long num, long prec)

.. function:: void zeta_ui_vec_odd(fmprb_ptr x, ulong start, long num, long prec)

    Computes `\zeta(s)` at num consecutive integers (respectively num
    even or num odd integers) beginning with `s = \mathrm{start} \ge 2`,
    automatically choosing an appropriate algorithm.


Euler-Maclaurin summation
-------------------------------------------------------------------------------

.. function:: void zeta_series_em_sum(fmpcb_ptr z, const fmpcb_t s, const fmpcb_t a, int deflate, ulong N, ulong M, long d, long prec)

.. function:: void zeta_series(fmpcb_ptr z, const fmpcb_t s, const fmpcb_t a, int deflate, long d, long prec)

    Evaluates the truncated Euler-Maclaurin sum of order `N, M` for the
    length-*d* truncated Taylor series of the Hurwitz zeta function
    `\zeta(s,a)` at `s`, using a working precision of *prec* bits.
    With `a = 1`, this gives the usual Riemann zeta function.

    If *deflate* is nonzero, `\zeta(s,a) - 1/(s-1)` is evaluated
    (which permits series expansion at `s = 1`).

    The *fmpcb_zeta_series* function chooses default values for `N, M`
    using *fmpcb_zeta_series_em_choose_param*,
    targeting an absolute truncation error of `2^{-\operatorname{prec}}`.

    The Euler-Maclaurin (EM) formula states that

    .. math ::

        \sum_{k=N}^U f(k) = \int_N^U f(t) dt + \frac{1}{2} \left(f(N) + f(U)\right)

                           + \sum_{k=1}^{M} \frac{B_{2k}}{(2k)!} \left( f^{(2k-1)}(U) - f^{(2k-1)}(N) \right)

                          - \int_N^U \frac{\tilde B_{2M}(t)}{(2M)!} f^{(2M)}(t) dt

    where `f` is a sufficiently differentiable function (for example,
    analytic), `B_n` is a Bernoulli number, and
    `\tilde B_n(t) = B_n(t-\lfloor t\rfloor)` is a periodic Bernoulli
    polynomial. If `f` decreases sufficiently rapidly, the formula
    remains valid after letting `U \to \infty`.

    To evaluate the Hurwitz zeta function, we set `f(k) = (a + k)^{-s}`,
    giving `\zeta(s,a) = \sum_{k=0}^{N-1} f(k) + \sum_{k=N}^{\infty} f(k)`,
    where EM summation is applied to the right sum.
    By choosing `M` and `N` large enough, and taking the standard
    logarithm branch cut, the EM formula gives an analytic
    continuation of `\zeta(s,a)` to all `a, s \in \mathbb{C}`
    (except for poles at `s = 1` and
    `\mathrm{Re}(s) > 0, a = 0, -1, -2, \ldots`). In order to
    evaluate derivatives with respect to `s` of `\zeta(s,a)`, we
    substitute `s \to s + x \in \mathbb{C}[[x]]`.

    We choose `N` such that `\Re(a+N) > 0`. Then the first integral is
    well-defined for `s` with `\Re(s) > 1` and has the closed form

    .. math ::

        \int_N^{\infty} f(t) dt = \int_N^{\infty}
            (a + t)^{-s}dt = \frac{(a+N)^{1-s}}{s-1},

    providing analytic continuation of this term with respect to `s`.
    Removing the singularity from this term also conveniently allows us
    to evaluate derivatives of `\zeta(s,a) - 1/(s-1)` at `s = 1`.

    The derivatives of `f(k)` are given by

    .. math ::

        f^{(r)}(k) = \frac{(-1)^r (s)_{r}}{(a+k)^{s+r}}

    where `(s)_{r} = s (s+1) \cdots (s+r-1)` denotes a rising factorial.
    Thus, the remainder integral becomes

    .. math ::

        R(s) = \int_N^{\infty} \frac{\tilde B_{2M}(t)}{(2M)!} \frac{(s)_{2M}}{(a+t)^{s+2M}} dt,

    valid when `\Re(a+N) > 0` and `\Re(s+2M-1) > 0`. We will use the
    stronger condition `\Re(a+N) > 1`.

    If `F = \sum_k f_k x^k \in \mathbb{C}[[x]]`, define
    `|F| = \sum_k |f_k| x^k` and `|F| \le |G|` if
    `\forall_k : |f_k| \le |g_k|`. It is easy to check that
    `|F + G| \le |F| + |G|` and `|FG| \le |F||G|`. With this notation,

    .. math ::

        |R(s+x)| = \left|\int_N^{\infty} \frac{\tilde B_{2M}(t)}{(2M)!}
            \frac{(s+x)_{2M}}{(a+t)^{s+x+2M}} dt\right|

        \le \int_N^{\infty} \left| \frac{\tilde B_{2M}(t)}{(2M)!}
            \frac{(s+x)_{2M}}{(a+t)^{s+x+2M}} \right| dt

        \le \frac{4 \left| (s+x)_{2M} \right|}{(2 \pi)^{2M}}
            \int_N^{\infty} \left| \frac{1}{(a+t)^{s+x+2M}} \right| dt

    where the fact that `|\tilde B_{2M}(x)| < 4 (2M)! / (2\pi)^{2M}` has
    been invoked. Thus it remains to bound the coefficients `R_k` satisfying

    .. math ::

        \int_N^{\infty} \left| \frac{1}{(a+t)^{s+x+2M}} \right| dt = 
            \sum_k R_k x^k, \quad
            R_k = \int_N^{\infty} \frac{1}{k!}
            \left| \frac{\log(a+t)^k}{(a+t)^{s+2M}} \right| dt.

    Writing `a = \alpha + \beta i`, where by assumption
    `\alpha + t \ge \alpha + N \ge 1`, we have

    .. math ::

        |\log(\alpha + \beta i + t)|
            = \left|\log(\alpha + t) +
            \log\left( 1 + \frac{\beta i}{\alpha + t}\right) \right|
            \le \log(\alpha + t) + \left|\log\left(1 + \frac{\beta i}{\alpha+t}\right)\right|

        = \log(\alpha+t) + \left|\frac{1}{2}\log\left(1+\frac{\beta^2}{(\alpha+t)^2}\right)
        + i\tan^{-1}\left(\frac{\beta}{\alpha + t}\right)\right| \le \log(\alpha + t) + C

    where

    .. math ::

        C = \frac{1}{2}\log\left(1+\frac{\beta^2}{(\alpha+N)^2}\right) +
            \tan^{-1}\left(\frac{|\beta|}{\alpha+N}\right) \le \frac{\beta^2}{2 (\alpha+N)^2}
            + \frac{|\beta|}{(\alpha+N)}.

    Also writing `s = \sigma + \tau i`, where by assumption `\sigma + 2M > 1`,
    we have

    .. math ::

        \frac{1}{|(\alpha+\beta i+t)^{\sigma+\tau i + 2M}|}
        = \frac{e^{\tau \operatorname{arg}(\alpha+\beta i+t)}}{|\alpha+\beta i+t|^{\sigma+2M}}
        \le \frac{K}{(\alpha + t)^{\sigma+2M}}

    where `K = \exp(\max(0, \tau \tan^{-1}(\beta / (\alpha + N))))`. Finally,

    .. math ::

        R_k \le \frac{K}{k!} \, I_k(N+\alpha, \sigma + 2M, C)

    with the `K` and `C` defined above, where `I_k(A,B,C)` denotes the
    sequence of integrals

    .. math ::

        I_k(A,B,C) \equiv \int_A^{\infty} t^{-B} (C + \log t)^k dt

    which can be evaluated as

    .. math ::

        I_k(A,B,C) = \frac{L_k}{(B-1)^{k+1} A^{B-1}}

    where `L_0 = 1`, `L_k = k L_{k-1} + D^k` and `D = (B-1) (C + \log A)`.

.. function:: void zeta_series_em_choose_param(fmpr_t bound, ulong * N, ulong * M, const fmpcb_t s, const fmpcb_t a, long d, long target, long prec)

    Chooses *N* and *M* using a default algorithm.


Power sums
-------------------------------------------------------------------------------

.. function:: void zeta_log_ui_from_prev(fmprb_t log_k1, ulong k1, fmprb_t log_k0, ulong k0, long prec)

    Computes `\log(k_1)`, given `\log(k_0)` where `k_0 < k_1`.
    At high precision, this function uses the formula
    `\log(k_1) = \log(k_0) + 2 \operatorname{atanh}((k_1-k_0)/(k_1+k_0))`,
    evaluating the inverse hyperbolic tangent using binary splitting
    (for best efficiency, `k_0` should be large and `k_1 - k_0` should
    be small). Otherwise, it ignores `\log(k_0)` and evaluates the logarithm
    the usual way.

.. function:: void zeta_powsum_series_naive(fmpcb_ptr z, const fmpcb_t s, const fmpcb_t a, long n, long len, long prec)

.. function:: void zeta_powsum_series_naive_threaded(fmpcb_ptr z, const fmpcb_t s, const fmpcb_t a, long n, long len, long prec)


    Computes

    .. math ::

        z = S(s,a,n) = \sum_{k=0}^{n-1} \frac{1}{(k+a)^{s+t}}

    as a power series in `t` truncated to length *len*. This function
    evaluates the sum naively term by term.
    The *threaded* version splits the computation
    over the number of threads returned by *flint_get_num_threads()*.

.. function:: void zeta_powsum_one_series_sieved(fmpcb_ptr z, const fmpcb_t s, long n, long len, long prec)

    Computes

    .. math ::

        z = S(s,1,n) \sum_{k=1}^n \frac{1}{k^{s+t}}

    as a power series in `t` truncated to length *len*.
    This function stores a table of powers that have already been calculated,
    computing `(ij)^r` as `i^r j^r` whenever `k = ij` is
    composite. As a further optimization, it groups all even `k` and
    evaluates the sum as a polynomial in `2^{-(s+t)}`.
    This scheme requires about `n / \log n` powers, `n / 2` multiplications,
    and temporary storage of `n / 6` power series. Due to the extra
    power series multiplications, it is only faster than the naive
    algorithm when *len* is small.

