.. _dlog:

**dlog.h** -- discrete logarithms mod ulong primes
===============================================================================

This module implements discrete logarithms, with the application
to Dirichlet characters in mind.

In particular, this module defines a :type:`dlog_precomp_t` structure
permitting to describe a discrete log problem  in some subgroup
of `\mathbb Z/p \mathbb Z` and store precomputed data for
faster computation of several such discrete logarithms.

When initializing this data, the user provides both a group description and the expected
number of subsequent discrete logarithms calls. The choice of algorithm and
the amount of stored data depend both on the structure of the group and this number.

No particular effort has been made towards single discrete logarithm
computation. Currently only machine size primepower moduli
are supported.

Types, macros and constants
-------------------------------------------------------------------------------

.. macro:: DLOG_NONE

   Return value when the discrete logarithm does not exist

.. type:: dlog_precomp_struct

.. type:: dlog_precomp_t

   Structure for discrete logarithm precomputed data.

.. function:: void dlog_precomp_clear(dlog_precomp_t pre)

Single evaluation
-------------------------------------------------------------------------------

.. function:: ulong dlog_once(ulong b, ulong a, const nmod_t mod, ulong n)

   Returns `x` such that `b = a^x` in `(\mathbb Z/mod \mathbb Z)^\times`,
   a has order *n*.

Precomputations
-------------------------------------------------------------------------------

.. function:: void dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)

   Precompute data for *num* discrete logarithms evaluations in the subgroup generated
   by *a* modulo *mod*, where *a* is known to have order *n*.

Specialized versions are available when specific information is known about the
group:

.. function:: void dlog_precomp_modpe_init(dlog_precomp_t pre, ulong a, ulong p, ulong e, ulong pe, ulong num)

   Assume that *a* generates the group of residues modulo `p^e`.

.. function:: void dlog_precomp_p_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong num)

   Assume that *a* has prime order *p*.

.. function:: void dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num)

   Assume that *a* has primepower order *p*.


Evaluation
-------------------------------------------------------------------------------

.. function:: ulong dlog_precomp(const dlog_precomp_t pre, ulong b)

   Returns `\log(b)` for the group described in *pre*

Vector evaluations
-------------------------------------------------------------------------------

.. function:: void dlog_vec_fill(ulong * v, ulong nv, ulong x)

   Sets values *v[k]* to *x* for all *k* less than *nv*.

.. function:: void dlog_vec_set_not_found(ulong * v, ulong nv, nmod_t mod)

   Sets values *v[k]* to :macro:`DLOG_NONE` for all *k* not coprime to *mod*.

.. function:: void dlog_vec(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   Sets *v[k]* to `\log(k,a)` times value *va*  for `0\leq k < nv`, where *a*
   has order *na*. *va* should be *1* for usual log computation.

.. function:: void dlog_vec_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)

   same parameters as before, but adds `\log(k,a)\times v_a`
   to *v[k]* and reduce modulo *order* instead of replacing the value. Indices
   *k* such that *v[k]* equals *DLOG_NONE* are ignored.

Algorithms
-------------------------------------------------------------------------------

Several discrete logarithms strategies are implemented:

- complete lookup table for small groups

- baby-step giant-step table

combined with mathematical reductions

- Pohlig-Hellman decomposition (chinese remainder decomposition on the
  order of the group and base `p` decomposition for primepower order)

- p-adic log for primepower modulus `p^e`.

For *dlog_vec* functions which compute the vector of discrete logarithms
of successive integers `1\dots n`:

- a simple loop on group elements avoiding all logarithms is done when
  the group size is comparable with the number of elements requested

- otherwise the logarithms are computed on primes and propagated by
  Eratosthene-like sieving on composite numbers.

- when several logarithms are already computed, a basic smoothing technique
  inspired by index-calculus is adopted to obtain larger logs from
  smaller ones.

- in the the present implementation, the full index-calculus method is not
  implemented.
