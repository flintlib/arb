.. _algorithms_hypergeometric:

Algorithms for hypergeometric functions
===============================================================================

.. _algorithms_hypergeometric_asymptotic_confluent:

Asymptotic series for the confluent hypergeometric function
-------------------------------------------------------------------------------

Let `U(a,b,z)` denote the confluent hypergeometric function of the second
kind with the principal branch cut, and
let `U^{*} = z^a U(a,b,z)`.
For all `z \ne 0` and `b \notin \mathbb{Z}` (but valid for all `b` as a limit),
we have (DLMF 13.2.42)

.. math ::

    U(a,b,z)
        = \frac{\Gamma(1-b)}{\Gamma(a-b+1)} M(a,b,z)
        + \frac{\Gamma(b-1)}{\Gamma(a)} z^{1-b} M(a-b+1,2-b,z).

Moreover, for all `z \ne 0` we have

.. math ::

    \frac{{}_1F_1(a,b,z)}{\Gamma(b)}
        = \frac{(-z)^{-a}}{\Gamma(b-a)} U^{*}(a,b,z)
        + \frac{z^{a-b} e^z}{\Gamma(a)} U^{*}(b-a,b,-z)

which is equivalent to DLMF 13.2.41 (but simpler in form).

We have the asymptotic expansion

.. math ::

    U^{*}(a,b,z) \sim {}_2F_0(a, a-b+1, -1/z)

where `{}_2F_0(a,b,z)` denotes a formal hypergeometric series, i.e.

.. math ::

    U^{*}(a,b,z) = \sum_{k=0}^{n-1} \frac{(a)_k (a-b+1)_k}{k! (-z)^k} + \varepsilon_n(z).

The error term `\varepsilon_n(z)` is bounded according to DLMF 13.7.
A case distinction is made depending on whether `z` lies in one
of three regions which we index by `R`.
Our formula for the error bound increases with the value of `R`, so we
can always choose the larger out of two indices if `z` lies in
the union of two regions.

Let `r = |b-2a|`.
If `\operatorname{Re}(z) \ge r`, set `R = 1`.
Otherwise, if `\operatorname{Im}(z) \ge r` or `\operatorname{Re}(z) \ge 0 \land |z| \ge r`, set `R = 2`.
Otherwise, if `|z| \ge 2r`, set `R = 3`.
Otherwise, the bound is infinite.
If the bound is finite, we have

.. math ::

    |\varepsilon_n(z)| \le 2 \alpha C_n \left|\frac{(a)_n (a-b+1)_n}{n! z^n} \right| \exp(2 \alpha \rho C_1 / |z|)

in terms of the following auxiliary quantities

.. math ::

    \sigma = |(b-2a)/z|

.. math ::

    C_n = \begin{cases}
    1                              & \text{if } R = 1 \\
    \chi(n)                        & \text{if } R = 2 \\
    (\chi(n) + \sigma \nu^2 n) \nu^n & \text{if } R = 3
    \end{cases}

.. math ::

    \nu = \left(\tfrac{1}{2} + \tfrac{1}{2}\sqrt{1-4\sigma^2}\right)^{-1/2} \le 1 + 2 \sigma^2

.. math ::

    \chi(n) = \sqrt{\pi} \Gamma(\tfrac{1}{2}n+1) / \Gamma(\tfrac{1}{2} n + \tfrac{1}{2})

.. math ::

    \sigma' = \begin{cases}
    \sigma & \text{if } R \ne 3 \\
    \nu \sigma & \text{if } R = 3
    \end{cases}

.. math ::

    \alpha = (1 - \sigma')^{-1}

.. math ::

    \rho = \tfrac{1}{2} |2a^2-2ab+b| + \sigma' (1+ \tfrac{1}{4} \sigma') (1-\sigma')^{-2}

.. _algorithms_hypergeometric_asymptotic_airy:

Asymptotic series for Airy functions
-------------------------------------------------------------------------------

Error bounds are based on Olver (DLMF section 9.7).
For `\arg(z) < \pi` and `\zeta = (2/3) z^{3/2}`, we have

.. math ::

    \operatorname{Ai}(z) = \frac{e^{-\zeta}}{2 \sqrt{\pi} z^{1/4}} \left[S_n(\zeta) + R_n(z)\right], \quad
    \operatorname{Ai}'(z) = -\frac{z^{1/4} e^{-\zeta}}{2 \sqrt{\pi}} \left[(S'_n(\zeta) + R'_n(z)\right]

.. math ::

    S_n(\zeta) = \sum_{k=0}^{n-1} (-1)^k \frac{u(k)}{\zeta^k}, \quad
    S'_n(\zeta) = \sum_{k=0}^{n-1} (-1)^k \frac{v(k)}{\zeta^k}

.. math ::

    u(k) = \frac{(1/6)_k (5/6)_k}{2^k k!}, \quad
    v(k) = \frac{6k+1}{1-6k} u(k).

Assuming that *n* is positive, the error terms are bounded by

.. math ::

    |R_n(z)|  \le C |u(n)| |\zeta|^{-n}, \quad |R'_n(z)| \le C |v(n)| |\zeta|^{-n}

where

.. math ::

    C = \begin{cases}
        2 \exp(7 / (36 |\zeta|)) & |\arg(z)| \le \pi/3 \\
        2 \chi(n) \exp(7 \pi / (72 |\zeta|)) & \pi/3 \le |\arg(z)| \le 2\pi/3 \\
        4 \chi(n) \exp(7 \pi / (36 |\operatorname{re}(\zeta)|)) |\cos(\arg(\zeta))|^{-n} & 2\pi/3 \le |\arg(z)| < \pi.
        \end{cases}

For computing Bi when *z* is roughly in the positive half-plane, we use the
connection formulas

.. math ::

    \operatorname{Bi}(z) = -i (2 w^{+1} \operatorname{Ai}(z w^{-2}) - \operatorname{Ai}(z))

    \operatorname{Bi}(z) = +i (2 w^{-1} \operatorname{Ai}(z w^{+2}) - \operatorname{Ai}(z))

where `w = \exp(\pi i/3)`. Combining roots of unity gives

.. math ::

    \operatorname{Bi}(z) = \frac{1}{2 \sqrt{\pi} z^{1/4}} [2X + iY]

.. math ::

    \operatorname{Bi}(z) = \frac{1}{2 \sqrt{\pi} z^{1/4}} [2X - iY]

.. math ::

    X = \exp(+\zeta) [S_n(-\zeta) + R_n(z w^{\mp 2})], \quad Y = \exp(-\zeta) [S_n(\zeta) + R_n(z)]

where the upper formula is valid
for `-\pi/3 < \arg(z) < \pi` and the lower formula is valid for `-\pi < \arg(z) < \pi/3`.
We proceed analogously for the derivative of Bi.

In the negative half-plane, we use the connection formulas

.. math ::

    \operatorname{Ai}(z) = e^{+\pi i/3} \operatorname{Ai}(z_1)  +  e^{-\pi i/3} \operatorname{Ai}(z_2)

.. math ::

    \operatorname{Bi}(z) = e^{-\pi i/6} \operatorname{Ai}(z_1)  +  e^{+\pi i/6} \operatorname{Ai}(z_2)

where `z_1 = -z e^{+\pi i/3}`, `z_2 = -z e^{-\pi i/3}`.
Provided that `|\arg(-z)| < 2 \pi / 3`, we have
`|\arg(z_1)|, |\arg(z_2)| < \pi`, and thus the asymptotic expansion
for Ai can be used. As before, we collect roots of unity to obtain

.. math ::

    \operatorname{Ai}(z) = A_1 [S_n(i \zeta)  + R_n(z_1)]
                         + A_2 [S_n(-i \zeta) + R_n(z_2)]

.. math ::

    \operatorname{Bi}(z) = A_3 [S_n(i \zeta)  + R_n(z_1)]
                         + A_4 [S_n(-i \zeta) + R_n(z_2)]

where `\zeta = (2/3) (-z)^{3/2}` and

.. math ::

    A_1 = \frac{\exp(-i (\zeta - \pi/4))}{2 \sqrt{\pi} (-z)^{1/4}}, \quad
    A_2 = \frac{\exp(+i (\zeta - \pi/4))}{2 \sqrt{\pi} (-z)^{1/4}}, \quad
    A_3 = -i A_1, \quad
    A_4 = +i A_2.

The differentiated formulas are analogous.

