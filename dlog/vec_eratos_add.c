/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "dlog.h"

/* assume non invertible and 1 mod n already set */
void
dlog_vec_eratos_add(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong p, k, n;
    dlog_precomp_t pre;
    n_primes_t iter;

    /* discrete log on primes */
    n = (nv < mod.n) ? nv : mod.n;
    dlog_precomp_n_init(pre, a, mod.n, na, n_prime_pi(n));

    n_primes_init(iter);
    while ((p = n_primes_next(iter)) < n)
    { 
        ulong wp, pe; 
        if (v[p] == DLOG_NOT_FOUND)
            continue; /* won't be attained another time */
        wp = nmod_mul(dlog_precomp(pre, p), va, order);
        /* fixme: could be faster sieving m*pe? but cannot
         * use v[p*m]=v[p]*v[m]... */
        for (pe = p; pe < n; pe *= p)
            for (k = pe; k < n; k += pe)
                if (v[k] != DLOG_NOT_FOUND)
                    v[k] = nmod_add(v[k],  wp, order);
    }
    n_primes_clear(iter);
    for (k = mod.n + 1; k < nv; k++)
        v[k] = v[k - mod.n];
}
