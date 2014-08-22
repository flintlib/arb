.. _algorithmselem:

Algorithms for elementary functions
===============================================================================

(This section is incomplete.)

We typically compute elementary functions using the following steps:
reduction to a combination of
standard function (exp, log, atan, sin, cos of a real argument), reduction
to a standard domain, convergence-accelerating argument reduction
(using functional equations and possibly precomputed lookup tables),
followed by Taylor series evaluation.

Arctangents
-------------------------------------------------------------------------------

It is sufficient to consider `x \in [0,1)`, since
`\operatorname{atan}(x) = \operatorname{atan}(-x)` for `x < 0`,
`\operatorname{atan}(1) = \pi/4`, and
`\operatorname{atan}(x) = \pi/2 - \operatorname{atan}(1/x)` for `x > 1`.

For sufficiently small `x`, we use the Taylor series

.. math ::

    \operatorname{atan}(x) = x - \frac{x^3}{3} + \frac{x^5}{5} - \ldots

Applying the argument-halving formula

.. math ::

    \operatorname{atan}(x) = 2 \operatorname{atan}\left(\frac{x}{1+\sqrt{1+x^2}}\right)

`r` times gives `x \le 2^{-r}`.

Applying the formula

.. math ::

    \operatorname{atan}(x) = \operatorname{atan}(p/q) +
        \operatorname{atan}(w),
        \quad w = \frac{qx-p}{px+q},
        \quad p = \lfloor qx \rfloor

gives `0 \le w < 1/q`. At low precision, picking a moderately large `q` (say `q = 2^8`),
and using a lookup table for
`\operatorname{atan}(p/q), p = 0, 1, \ldots, q-1`, is much better than repeated argument-halving.
This transformation can be applied repeatedly with a sequence of increasing
values of `q`. For example, `(q_1 = 2^4, q_2 = 2^8)` requires
only `2^5` precomputed table entries, but the evaluation
costs one extra division).

At high precision, the `\operatorname{atan}(p/q)` values
can be evaluated using binary splitting. Choosing `q = 2, 4, 8, \ldots`
results in a version of the bit-burst a algorithm.

Error propagation for arctangents
-------------------------------------------------------------------------------

A generic derivative-based error bound is

.. math ::

    \sup_{\xi \in [-1,1]} |\operatorname{atan}(m) - \operatorname{atan}(m+\xi r)|
        \le \frac{r}{1+\max(0,|m|-r)^2} \le r.

An exact representation for the propagated error is given by

.. math ::

    \sup_{\xi \in [-1,1]}  |\operatorname{atan}(m) - \operatorname{atan}(m + \xi r)| =
        \begin{cases}
        \operatorname{atan}\left(\frac{r}{1+|m|(|m|-r)}\right) & \text{if } r \le |m| \\
        \frac{\pi}{2} - \operatorname{atan}\left(\frac{1 + |m| (|m|-r)}{r}\right) & \text{if } r > |m|
        \end{cases}.

Logarithms
-------------------------------------------------------------------------------

It is sufficient to consider `\log(1+x)` where `x \in [0,1)`, since
`\log(t 2^n) = \log(t) + n \log(2)`.

We only use the Taylor series

.. math ::

    \log(1+x) = x - \frac{x^2}{2} + \frac{x^3}{3} - \ldots

directly when `x` is so small that the linear or quadratic
term gives full accuracy. In general we use the more efficient series

.. math ::

    \operatorname{atanh}(x) = x + \frac{x^3}{3} + \frac{x^5}{5} + \ldots

together with the identity `\log(1+x) = 2 \operatorname{atanh}(x/(2+x))`.

Applying the argument-halving formula

.. math ::

    \log(1+x) = 2 \log\left(\sqrt{1+x}\right)

`r` times gives `x \le 2^{-r}`.

Applying the formula

.. math ::

    \log(1+x) = \log(p/q) +
        \log(1+w),
        \quad w = \frac{qx-p}{p+q},
        \quad p = \lfloor qx \rfloor

gives `0 \le w < 1/q` (see analogous remarks for the arctangent).

Error propagation for logarithms
-------------------------------------------------------------------------------

A generic derivative-based error bound is

.. math ::

    \sup_{\xi \in [-1,1]} |\log(m) - \log(m+\xi r)|
        \le \frac{r}{m-r}.

An exact representation for the propagated error is given by

.. math ::

    \sup_{\xi \in [-1,1]} |\log(m) - \log(m+\xi r)| = \log(1 + r/(m-r)).

Of course, these formulas require `m > r \ge 0` (otherwise, the real
logarithm is undefined).

