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
#include <math.h>

#define vbs 1
#define FACTOR_RATIO 4

static int
factor_until(ulong * n, ulong nlim, const ulong * p, ulong pmax, ulong * fp, int * fe)
{
    int i, j;
    for (i = 0, j = 0; *n >= nlim && p[j] < pmax; j++)
    {
        int e = n_remove(n, p[j]);
        if (e)
        {
            fp[i] = p[j];
            fe[i] = e;
            i++;
        }
    }
    return i;
}

ulong
dlog_vec_pindex_factorgcd(ulong * v, ulong nv, ulong p, nmod_t mod, ulong a, ulong na, ulong loga, ulong logm1, nmod_t order, int maxtry)
{
    int nm = 0, ng = 0;
    ulong pm, logm, pmax;
    ulong u[2], r[2], t;
    ulong up[15], rp[15];
    int ue[15], re[15];
    const ulong * prime;
    prime = n_primes_arr_readonly(p);
    pmax = p / FACTOR_RATIO;
    pm = p;
    logm = 0;
    while (nm++ < maxtry)
    {
        int i, j, iu, ir;
        ulong logr;
        pm = nmod_mul(pm, a, mod);
        logm = nmod_add(logm, loga, order);
        /*
           if (2 * pm > mod.n)
           {
           pm = nmod_neg(pm, mod);
           logm = nmod_add(logm, logm1, order);
           }
           */
        /* half gcd u * pm + v * mod = r, ignore v */
        u[0] = 0; r[0] = mod.n;
        u[1] = 1; r[1] = pm;
        i = 1; j = 0; /* flip flap */
        while (r[i] > u[i])
        {
            ng++;
            if (r[i] < nv && v[r[i]] != NOT_FOUND && u[i] < nv && v[u[i]] != NOT_FOUND)
            {
                /* early smooth detection: occurs for primes < 30 bits */
                ulong x;
                /* chi(-1)^j*chi(u)*chi(p)*chi(m)=chi(r) */
                x = nmod_sub(v[r[i]], nmod_add(v[u[i]], logm, order), order);
                if (j)
                    x = nmod_add(x, logm1, order);
                flint_printf("[sieve early] %wu * %wu^%wu = %wu [%wu]\n",
                        p, a, nm, pm, mod.n);
                flint_printf("[sieve early] found %wu * %wu = (-1)^%d*%wu [%wu]\n",
                        u[i], pm, j, r[i], mod.n);
                flint_printf("[sieve early] log(%wu^%wu) = %wu * %wu = %wu [%wu]\n",
                        a, nm, nm, loga, logm, order.n);
                flint_printf("[    on logs] %wu + %wu + log(%wu) = %d * %wu + %wu [%wu]\n",
                        v[u[i]],logm,p,j,logm1,v[r[i]], order.n);
                flint_printf("[    hence  ] log(%wu) = %wu\n", p, x);
                return x;
            }
            j = i; i = 1 - i; /* switch */
            t = r[i] / r[j];
            r[i] = r[i] % r[j];
            u[i] = u[i] + t * u[j]; /* times (-1)^j */
        };
        flint_printf("[sieve factor] %wu * %wu^%wu = %wu [%wu]\n",
                p, a, nm, pm, mod.n);
        flint_printf("[sieve factor] found %wu * %wu = (-1)^%d*%wu [%wu]\n",
                u[i], pm, j, r[i], mod.n);
        flint_printf("[sieve factor] logm = %wu [A=%wu,logA=%wu,nm=%wu]\n",
                logm,a,loga,nm);
        logr = (j) ? logm1 : 0;
        flint_printf("[sieve factor] logr = %wu\n",logr);
        /* try to factor both r[i] and u[i] */
        iu = factor_until(&u[i], nv, prime, pmax, up, ue);
        if (u[i] >= nv || v[u[i]] == NOT_FOUND)
            continue;
        flint_printf("[sieve factor] u: found %d factors up to %wu\n",iu,u[i]);
        ir = factor_until(&r[i], nv, prime, pmax, rp, re);
        if (r[i] >= nv || v[r[i]] == NOT_FOUND)
            continue;
        flint_printf("[sieve factor] r: found %d factors up to %wu\n",ir,r[i]);
        /* log(u)+log(p)+log(m)=log(r) */
        logm = nmod_add(logm, v[u[i]], order);
        logr = nmod_add(logr, v[r[i]], order);
        for (i=0; i < ir; i++)
            logr = nmod_add(logr, (re[i] * v[rp[i]]) % order.n, order);
        flint_printf("[sieve factor] logr = %wu\n",logr);
        for (i=0; i < iu; i++)
            logm = nmod_add(logm, (ue[i] * v[up[i]]) % order.n, order);
        flint_printf("[sieve factor] logm = %wu\n",logm);
        flint_printf("[sieve factor] log(%wu^%wu) = %wu * %wu = %wu [%wu]\n",
                a, nm, nm, loga, logm, order.n);
        flint_printf("[    on logs] %wu + log(%wu) = %d * %wu + %wu [%wu]\n",
                logm,p,j,logm1, logr, order.n);
        flint_printf("[    hence  ] log(%wu) = %wu\n", p, nmod_sub(logr, logm, order));

        return nmod_sub(logr, logm, order);
    }
    return NOT_FOUND;
}
