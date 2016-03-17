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

void
dlog_vec_eratos_subgroup(ulong *v, ulong nv, ulong a, ulong va, ulong M, nmod_t mod, ulong na, nmod_t order)
{
    ulong p, pmax;
    dlog_precomp_t pre;
    n_primes_t iter;

    /* discrete log on primes */
    pmax = (nv < mod.n) ? nv : mod.n;
    dlog_precomp_n_init(pre, a, mod.n, na, n_prime_pi(nv));

    n_primes_init(iter);
    while ((p = n_primes_next(iter)) < pmax)
    { 
        ulong k, pM, wp; 
        if (mod.n % p == 0)
            continue; /* won't be attained another time */
        pM = (M) ? nmod_pow_ui(p, M, mod) : p;
        wp = nmod_mul(dlog_precomp(pre, pM), va, order);
        for (pM = p; pM < nv; pM *= p)
            for (k = pM; k < nv; k += pM)
                v[k] = nmod_add(v[k],  wp, order);
    }
    n_primes_clear(iter);
}
