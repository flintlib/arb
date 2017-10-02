.. _acb-dft:

**acb_dft.h** -- Discrete Fourier Transform on finite abelian groups
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

Let *G* be a finite abelian group, and `\chi` a character of *G*.
For any map `f:G\to\mathbb C`, the discrete fourier transform
`\hat f:\hat G\to \mathbb C` is defined by

.. math::

   \hat f(\chi) = \sum_{x\in G}\overline{\chi(x)}f(x)


Fast Fourier Transform techniques allow to compute efficiently
all values `\hat f(\chi)` by reusing common computations.

Specifically, if `H\triangleleft G` is a subgroup of size `M` and index
`[G:H]=m`, then writing `f_x(h)=f(xh)` the translate of `f` by representatives
`x` of `G/H`, one has a decomposition

.. math::
   
   \hat f(\chi) = \sum_{x\in G/H} \overline{\chi(x)} \hat{f_x}(\chi_{H})

so that the DFT on `G` can be computed using `m` DFT  on `H` (of
appropriate translates of `f`), then `M` DFT on `G/H`, one for
each restriction `\chi_{H}`.

This decomposition can be done recursively.

Note that by inversion formula

.. math::

   \widehat{\hat f}(\chi) = \#G\times f(\chi^{-1})

it is straightforward to recover `f` from its DFT `\hat f`.

DFT on Z/nZ
-------------------------------------------------------------------------------

If `G=\mathbb Z/n\mathbb Z`, we compute the DFT according to the usual convention

.. math::

   w_x = \sum_{y\bmod n} v_y e^{-\frac{2iπ}nxy}

.. function:: void acb_dirichlet_dft_naive(acb_ptr w, acb_srcptr v, slong n, slong prec)

.. function:: void acb_dirichlet_dft_crt(acb_ptr w, acb_srcptr v, slong n, slong prec)

.. function:: void acb_dirichlet_dft_cyc(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Set *w* to the DFT of *v* of length *len*.

   The first variant uses the naive `O(n^2)` algorithm.
   
   The second one uses CRT to express `Z/nZ` as a product of cyclic groups.

   The *cyc* version uses each prime factor of `m` of `n` to decompose with
   the subgroup `H=m\mathbb Z/n\mathbb Z`.


DFT on products
-------------------------------------------------------------------------------

A finite abelian group is isomorphic to a product of cyclic components

.. math::

   G = \bigoplus_{i=1}^r \mathbb Z/n_i\mathbb Z

then a character is a product of characters of all components and the DFT reads

.. math::

   \hat f(x_1,\dots x_r) = \sum_{y_1\dots y_r} f(y_1,\dots y_r)
   e^{-2iπ\sum\frac{x_i y_i}{n_i}}

We assume that `f` is given by a vector of length `\prod n_i` corresponding
to a lexicographic ordering of the values `y_1,\dots y_r`, and the computation
returns the same indexing for values of `\hat f`.

.. function:: void acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec);

   Computes the DFT on the group product of *num* cyclic components of sizes *cyc*.

Precomputations
-------------------------------------------------------------------------------

Convolution
-------------------------------------------------------------------------------

For functions `f` and `g` on `G` we consider the convolution

.. math::

   (f \star g)(x) = \sum_{y\in G} f(x-y)g(y)

which satisfies

.. math::

   \widehat{f \star g}(\chi) = \hat f(\chi)\hat g(\chi)


