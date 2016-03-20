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

#define SIEVE_START 100

/* assume v[k] = -1 for bad primes? */
/* loop on small components and if needed keep one subgroup for DLOG + sieve */
  void
dlog_vec_crt(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    n_factor_t fac;
    ulong maxloop;
    int k;

    dlog_vec_fill(v, nv, 0);
    maxloop = LOOP_MAX_FACTOR * nv;
    n_factor_init(&fac);
    n_factor(&fac, na, 1);
    dlog_n_factor_group(&fac, maxloop);
    for (k = 0; k < fac.num; k++)
    {
        ulong m, M, aM, uM, vaM;
        m = fac.p[k];
        M = na / m;
        aM = nmod_pow_ui(a, M, mod);
        uM = M * n_invmod(M % m, m); /* uM < n */
        vaM = nmod_mul(va, uM % order.n, order);
        if (fac.exp[k] == G_SMALL)
            dlog_vec_loop_subgroup(v, nv, aM, vaM, mod, m, order);
        else
        {
            if (nv <= SIEVE_START)
                dlog_vec_eratos_subgroup(v, nv, aM, vaM, M, mod, m, order);
            else
                dlog_vec_sieve_subgroup(v, nv, aM, vaM, M, mod, m, order);
        }
    }
}
