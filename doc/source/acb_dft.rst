.. _acb-dft:

**acb_dft.h** -- Discrete Fourier Transform on finite abelian groups
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

Discrete Fourier Transform
-------------------------------------------------------------------------------

Let *G* be a finite abelian group, and `\chi` a character of *G*.
For any map `f:G\to\mathbb C`, the discrete fourier transform
`\hat f:\hat G\to \mathbb C` is defined by

.. math::

   \hat f(\chi) = \sum_{x\in G}\overline{\chi(x)}f(x)

Note that by inversion formula

.. math::

   \widehat{\hat f}(\chi) = \#G\times f(\chi^{-1})

it is straightforward to recover `f` from its DFT `\hat f`.

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

DFT on Z/nZ
-------------------------------------------------------------------------------

If `G=\mathbb Z/n\mathbb Z`, we compute the DFT according to the usual convention

.. math::

   w_x = \sum_{y\bmod n} v_y e^{-\frac{2iπ}nxy}

.. function:: void acb_dft_naive(acb_ptr w, acb_srcptr v, slong n, slong prec)

.. function:: void acb_dft_crt(acb_ptr w, acb_srcptr v, slong n, slong prec)

.. function:: void acb_dft_cyc(acb_ptr w, acb_srcptr v, slong n, slong prec)

.. function:: void acb_dft_bluestein(acb_ptr w, acb_srcptr v, slong n, slong prec)

.. function:: void acb_dft(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Set *w* to the DFT of *v* of length *len*.

   The first variant uses the naive `O(n^2)` algorithm.

   The second one uses CRT to express `Z/nZ` as a product of cyclic groups.

   The *cyc* version uses each prime factor of `m` of `n` to decompose with
   the subgroup `H=m\mathbb Z/n\mathbb Z`.

   The *bluestein* version converts the computation to a radix 2 one using Bluestein's convolution trick.

   The default version uses an automatic choice of algorithm (in most cases *crt*).

.. function:: void acb_dft_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec)

   Compute the inverse DFT of *v* into *w*.

.. function:: void acb_dft_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)

.. function:: void acb_dft_rad2_inplace(acb_ptr w, int e, slong prec)

   Computes the DFT of *v* into *w*, where *v* and *w* have size $2^e$,
   using a radix 2 FFT. The ``inplace`` version corresponds to *v=w*.

.. function:: void acb_dft_inverse_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)

.. function:: void acb_dft_inverse_rad2_inplace(acb_ptr w, int e, slong prec)

   Computes the inverse DFT of *v* into *w*, where *v* and *w* have size $2^e$,
   using a radix 2 FFT. The ``inplace`` version corresponds to *v=w*.

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

.. function:: void acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec)

   Computes the DFT on the group product of *num* cyclic components of sizes *cyc*.

Precomputations
-------------------------------------------------------------------------------

If several computations are to be done on the same group, the FFT scheme
should be reused.

.. type:: acb_dft_pre_struct

.. type:: acb_dft_pre_t

    Stores a fast DFT scheme on :math:`\mathbb Z/n\mathbb Z`
    as a recursive decomposition into simpler DFT
    with some tables of roots of unity.

    An *acb_dft_pre_t* is defined as an array of *acb_dft_pre_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec)

   Initializes the fast DFT scheme of length *len*, using an automatic choice of
   algorithms depending on the factorization of *len*.

   The length *len* is stored as *pre->n*.

.. function:: void acb_dft_precomp_clear(acb_dft_pre_t pre)

   Clears *pre*.

.. function:: void acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)

   Computes the DFT of the sequence *v* into *w* by applying the precomputed scheme
   *pre*. Both *v* and *w* must have length *pre->n*.

.. function:: void acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)

   Compute the inverse DFT of *v* into *w*.

Convolution
-------------------------------------------------------------------------------

For functions `f` and `g` on `G` we consider the convolution

.. math::

   (f \star g)(x) = \sum_{y\in G} f(x-y)g(y)

.. function:: void acb_dft_convol_naive(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)

.. function:: void acb_dft_convol_rad2(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)

.. function:: void acb_dft_convol(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)

   Sets *w* to the convolution of *f* and *g* of length *len*.

   The *naive* version simply uses the definition.

   The *rad2* version embeds the sequence into a power of 2 length and
   uses the formula

   .. math::

      \widehat{f \star g}(\chi) = \hat f(\chi)\hat g(\chi)

   to compute it using three radix 2 FFT.

   The default version uses radix 2 FFT unless *len* is a product of small
   primes where a non padded FFT is faster.
