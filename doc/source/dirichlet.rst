.. _dirichlet:

**dirichlet.h** -- Dirichlet characters
===================================================================================

*Warning: the interfaces in this module are experimental and may change
without notice.*

This module allows working with Dirichlet characters algebraically.
For evaluations of characters as complex numbers and Dirichlet L-functions, 
see the :ref:`acb_dirichlet` module.

Multiplicative group modulo *q*
-------------------------------------------------------------------------------

Working with Dirichlet characters mod *q* consists mainly
in going from residue classes mod *q* to exponents on a set
of generators of the group.

This implementation relies on the Conrey numbering scheme
introduced in the LMFDB, which is an explicit choice of isomorphism

.. math::

   (\mathbb Z/q\mathbb Z)^\times & \to &\bigoplus_i \mathbb Z/\phi_i\mathbb Z \\
   x & \mapsto & (e_i)

We call *number* a residue class `x` modulo *q*, and *log* the
corresponding vector `(e_i)` of exponents of Conrey generators.

Going from a *log* to the corresponding *number* is a cheap
operation called exp, while the converse requires computing discrete
logarithms.

.. type:: dirichlet_group_struct

.. type:: dirichlet_group_t

    Represents the group of Dirichlet characters mod *q*.

    An *dirichlet_group_t* is defined as an array of *dirichlet_group_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void dirichlet_group_init(dirichlet_group_t G, ulong q)

    Initializes *G* to the group of Dirichlet characters mod *q*.

    This method computes a canonical decomposition of *G* in terms of cyclic
    groups, which are the mod `p^e` subgroups for `p^e\|q`.
    In particular *G* contains:

    - the number *num* of components

    - the generators

    - the exponent *expo* of the group

    It does *not* automatically precompute lookup tables
    of discrete logarithms or numerical roots of unity, and can therefore
    safely be called even with large *q*.

    For implementation reasons, the largest prime factor of *q* must not
    exceed `10^{12}` (an abort will be raised). This restriction could
    be removed in the future.

.. function:: void dirichlet_subgroup_init(dirichlet_group_t H, const dirichlet_group_t G, ulong h)

   Given an already computed group *G* mod `q`, initialize its subgroup *H*
   defined mod `h\mid q`. Precomputed discrete log tables are inherited.

.. function:: void dirichlet_group_clear(dirichlet_group_t G)

    Clears *G*. Remark this function does *not* clear the discrete logarithm
    tables stored in *G* (which may be shared with another group).

.. function:: void dirichlet_group_dlog_precompute(dirichlet_group_t G, ulong num)

    Precompute decomposition and tables for discrete log computations in *G*,
    so as to minimize the complexity of *num* calls to discrete logarithms.

    If *num* gets very large, the entire group may be indexed.

.. function:: void dirichlet_group_dlog_clear(dirichlet_group_t G, ulong num)

   Clear discrete logarithm tables in *G*. When discrete logarithm tables are
   shared with subgroups, those subgroups must be cleared before clearing the
   tables.

Conrey elements
-------------------------------------------------------------------------------

.. type:: dirichlet_conrey_struct

.. type:: dirichlet_conrey_t

    Represents elements of the unit group mod *q*, keeping both the
    *number* (residue class) and *log* (exponents on the group
    generators).

.. function:: void dirichlet_conrey_log(dirichlet_conrey_t x, const dirichlet_group_t G, ulong m)

    Sets *x* to the element of number *m*, computing its log using discrete
    logarithm in *G*.

.. function:: ulong dirichlet_conrey_exp(dirichlet_conrey_t x, const dirichlet_group_t G)

    Compute the reverse operation.

.. function:: void dirichlet_conrey_one(dirichlet_conrey_t x, const dirichlet_group_t G)

    Sets *x* to the *number* `1\in G`, having *log* `[0,\dots 0]`.

.. function:: void dirichlet_conrey_first_primitive(dirichlet_conrey_t x, const dirichlet_group_t G)

    Sets *x* to the first primitive element of *G*, having *log* `[1,\dots 1]`,
    or `[0, 1, \dots 1]` if `8\mid q`.

.. function:: void dirichlet_conrey_set(dirichlet_conrey_t x, const dirichlet_group_t G, const dirichlet_conrey_t y)

    Sets *x* to the element *y*.

.. function:: int dirichlet_conrey_next(dirichlet_conrey_t x, const dirichlet_group_t G)

    Sets *x* to the next conrey element in *G* with lexicographic ordering.

    The return value
    is the index of the last updated exponent of *x*, or *-1* if the last
    element has been reached.

    This function allows to iterate on the elements of *G* looping on their *log*.
    Note that it produces elements in seemingly random *number* order.

    The following template can be used to loop over all elements *x* in *G*::

        acb_conrey_one(x, G);
        do {
            /* use Conrey element x */
        } while (dirichlet_conrey_next(x, G) >= 0);

.. function:: int dirichlet_conrey_next_primitive(dirichlet_conrey_t x, const dirichlet_group_t G)

    Same as :func:`dirichlet_conrey_next`, but jumps to the next element
    corresponding to a primitive character of *G*.

.. function:: ulong dirichlet_index_conrey(const dirichlet_group_t G, const dirichlet_conrey_t x);

    Returns the lexicographic index of *x* as an integer in `0\dots \varphi(q)`.

.. function:: void dirichlet_conrey_index(dirichlet_conrey_t x, const dirichlet_group_t G, ulong j)

    Sets *x* to the Conrey element of lexicographic index *j*.

.. function:: int dirichlet_conrey_eq(const dirichlet_conrey_t x, const dirichlet_conrey_t y)

.. function:: int dirichlet_conrey_eq_deep(const dirichlet_group_t G, const dirichlet_conrey_t x, const dirichlet_conrey_t y)

   Return 1 if *x* equals *y*.
   The second version checks every byte of the representation and is intended for testing only.

Dirichlet characters
-------------------------------------------------------------------------------

Dirichlet characters take value in a finite cyclic group of roots of unity plus zero.

When evaluation functions return a *ulong*, this number corresponds to the
power of a primitive root of unity, the special value *DIRICHLET_CHI_NULL*
encoding the zero value.

The Conrey numbering scheme makes explicit the mathematical fact that
the group *G* is isomorphic to its dual, so that a character is described by
a *number*.

.. math::

   \begin{array}{ccccc}
   (\mathbb Z/q\mathbb Z)^\times \times (\mathbb Z/q\mathbb Z)^\times & \to & \bigoplus_i \mathbb Z/\phi_i\mathbb Z \times \mathbb Z/\phi_i\mathbb Z & \to &\mathbb C \\
   (m,n) & \mapsto& (a_i,b_i) &\mapsto& \chi_q(m,n) = \exp(2i\pi\sum \frac{a_ib_i}{\phi_i} )
   \end{array}

.. function:: ulong dirichlet_ui_pairing(const dirichlet_group_t G, ulong m, ulong n)

.. function:: ulong dirichlet_ui_pairing_conrey(const dirichlet_group_t G, const dirichlet_conrey_t a, const dirichlet_conrey_t b)

   Compute the value of the Dirichlet pairing on numbers *m* and *n*, as
   exponent modulo *G->expo*.
   The second form takes the Conrey index *a* and *b*, and does not take discrete
   logarithms.

   The returned value is the numerator of the actual value exponent mod the group exponent *G->expo*.

Character type
-------------------------------------------------------------------------------

.. type:: dirichlet_char_struct

.. type:: dirichlet_char_t

    Represents a Dirichlet character. This structure contains various
    useful invariants such as the order, the parity and the conductor of the character.

    An *dirichlet_char_t* is defined as an array of *dirichlet_char_struct*
    of length 1, permitting it to be passed by reference.

.. function:: void dirichlet_char_init(dirichlet_char_t chi, const dirichlet_group_t G)

    Initializes *chi* to an element of the group *G* and sets its value
    to the principal character.

.. function:: void dirichlet_char_clear(dirichlet_char_t chi)

    Clears *chi*.

.. function:: void dirichlet_char(dirichlet_char_t chi, const dirichlet_group_t G, ulong n)

    Sets *chi* to the Dirichlet character of number *n*, using Conrey numbering scheme.
    This function performs a discrete logarithm in *G*.

.. function:: void dirichlet_char_conrey(dirichlet_char_t chi, const dirichlet_group_t G, const dirichlet_conrey_t x)

    Sets *chi* to the Dirichlet character corresponding to *x*.

.. function:: int dirichlet_char_eq(const dirichlet_char_t chi1, const dirichlet_char_t chi2)

.. function:: int dirichlet_char_eq_deep(const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2)

   Return 1 if *chi1* equals *chi2*.
   The second version checks every byte of the representation and is intended for testing only.

.. function:: int dirichlet_char_is_principal(const dirichlet_char_t chi)

    Return 1 if *chi* is the principal character mod *q*.

.. function:: void dirichlet_char_one(dirichlet_char_t chi, const dirichlet_group_t G)

    Sets *chi* to the principal character.

.. function:: void dirichlet_char_set(dirichlet_char_t chi1, const dirichlet_group_t G, const dirichlet_char_t chi2)

    Sets *chi1* to the character *chi2*.

.. function:: int dirichlet_char_next(dirichlet_char_t chi, const dirichlet_group_t G)

    Sets *x* to the next character in *G* with lexicographic Conrey ordering
    (see :func:`dirichlet_conrey_next`). The return value
    is the index of the last updated exponent of *x*, or *-1* if the last
    element has been reached.

.. function:: int dirichlet_char_next_primitive(dirichlet_char_t chi, const dirichlet_group_t G)

    Like :func:`dirichlet_char_next`, but only generates primitive
    characters.

Character properties
-------------------------------------------------------------------------------

As a consequence of the Conrey numbering, all these numbers are available at the
level of *number* and Conrey *log* elements, and for *char*.
No discrete log computation is performed.

.. function:: ulong dirichlet_number_primitive(const dirichlet_group_t G)

   Return the number of primitive elements in *G*.

.. function:: ulong dirichlet_ui_conductor(const dirichlet_group_t G, ulong a)

.. function:: ulong dirichlet_conrey_conductor(const dirichlet_group_t G, const dirichlet_conrey_t x)

.. function:: ulong dirichlet_char_conductor(const dirichlet_char_t chi)

   Return the *conductor* of `\chi_q(a,\cdot)`, that is the smallest `r` dividing `q`
   such `\chi_q(a,\cdot)` can be obtained as a character mod `r`.
   This number is precomputed for the *char* type.

.. function:: int dirichlet_ui_parity(const dirichlet_group_t G, ulong a)

.. function:: int dirichlet_conrey_parity(const dirichlet_group_t G, const dirichlet_conrey_t x)

.. function:: int dirichlet_char_parity(const dirichlet_char_t chi)

   Return the *parity* `\lambda` in `\{0,1\}` of `\chi_q(a,\cdot)`, such that
   `\chi_q(a,-1)=(-1)^\lambda`.
   This number is precomputed for the *char* type.

.. function:: ulong dirichlet_ui_order(const dirichlet_group_t G, ulong a)

.. function:: ulong dirichlet_conrey_order(const dirichlet_group_t G, const dirichlet_conrey_t x)

.. function:: ulong dirichlet_char_order(const dirichlet_char_t chi)

   Return the order of `\chi_q(a,\cdot)` which is the order of `a\bmod q`.
   This number is precomputed for the *char* type.

.. function:: int dirichlet_char_is_real(const dirichlet_char_t chi)

   Return 1 if *chi* is a real character (iff it has order `\leq 2`).

Character evaluation
-------------------------------------------------------------------------------

The image of a Dirichlet character is a finite cyclic group. Dirichlet
character evaluations are exponents in this group.

.. function:: ulong dirichlet_ui_chi_conrey(const dirichlet_group_t G, const dirichlet_char_t chi, const dirichlet_conrey_t x)

.. function:: ulong dirichlet_ui_chi(const dirichlet_group_t G, const dirichlet_char_t chi, ulong n)

   Compute that value `\chi(n)` as the exponent mod the order of `\chi`.

Vector evaluation
-------------------------------------------------------------------------------

.. function:: void dirichlet_ui_chi_vec(ulong * v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv)

   Compute the list of exponent values *v[k]* for `0\leq k < nv`.

Character operations
-------------------------------------------------------------------------------

.. function:: void dirichlet_conrey_mul(dirichlet_conrey_t c, const dirichlet_group_t G, const dirichlet_conrey_t a, const dirichlet_conrey_t b)

.. function:: void dirichlet_char_mul(dirichlet_char_t chi12, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2)

   Multiply two characters in the same group.

.. function:: void dirichlet_conrey_pow(dirichlet_conrey_t c, const dirichlet_group_t G, const dirichlet_conrey_t a, ulong n)

   Take the power of some character.

Implementation notes
-------------------------------------------------------------------------------

The current implementation introduces a *char* type which contains a *conrey*
log plus additional information which

- makes evaluation of a single character a bit faster

- has some initialization cost.

Even if it is straightforward to convert a *conrey* log to the
corresponding *char*, looping is faster at the
level of Conrey representation. Things can be improved on this aspect
but it makes code more intricate.
