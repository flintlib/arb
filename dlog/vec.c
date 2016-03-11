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

#include "math.h"
#include "dlog.h"
#define NOT_FOUND UWORD_MAX

/* vector of log(k,a)*loga % order in Z/modZ */ 
void
dlog_vec_loop(ulong * v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order)
{
    ulong x, xp;
    long vx = 0;
    for(x = a; x != 1; x = nmod_mul(x, a, mod))
    {
      vx = nmod_add(vx, va, order);
      for(xp = x; xp < nv; xp+=mod.n)
          v[xp] = nmod_add(v[xp], vx, order);
    }
}

static ulong
logp_sieve(ulong * v, ulong nv, ulong p, ulong mod, ulong logm1, const nmod_t order, int maxtry)
{
    int i, j, nr = 0, nm = 0, ng = 0;
    ulong l, pm, logm;
    ulong u[2], r[2], t;
#if rnd
    flint_rand_t state;
    flint_randinit(state);
#endif
#if vbs
    flint_printf("\nEnter logp_sieve p=%wu mod %wu...\n", p, mod);
#endif
    while (1) {
        /* find multiplier m */
        logm = 0;
        pm = l = p;
        do {
            nm++;
            /* random ? pb when p lies in a small subgroup */
            do {
                nr++;
#if rnd
                l = 1 + n_randint(state, p - 1);
#else
                l = (l > 1) ? l - 1 : p - 1;
#endif
            } while (v[l] == NOT_FOUND);
            pm *= l;
            logm += v[l];
        } while (pm < mod);
        pm = pm % mod;
#if vbs
        flint_printf("[pm=%wu, v[pm]=%ld]", pm, v[pm]);
#endif
        /* half gcd u * pm + v * mod = r, ignore v */
        u[0] = 0; r[0] = mod;
        u[1] = 1; r[1] = pm;
        i = 1; j = 0; /* flip flap */
        do {
            ng++;
#if vbs
            flint_printf("[r=%d, v[r]=%d, u=%d, v[u]=%d]\n",
                    r[i],v[r[i]], u[i], v[u[i]]);
#endif
            if (r[i] < nv && v[r[i]] != NOT_FOUND && u[i] < nv && v[u[i]] != NOT_FOUND)
            {
                ulong x;
                /* chi(-1)^j*chi(u)*chi(p)*chi(m)=chi(r) */
                x = nmod_sub(v[r[i]], nmod_add(v[u[i]], logm, order), order);
                if (j)
                    x = nmod_add(x, logm1, order);
#if rnd
                flint_randclear(state);
#endif
                return x;
            }

            j = i; i = 1 - i; /* switch */
            t = r[i] / r[j];
            r[i] = r[i] % r[j];
            u[i] = u[i] + t * u[j]; /* (-1)^j */

        } while (r[i] > 0 && u[i] < p);
        if (nm > maxtry)
            return NOT_FOUND;
    }
}

void
dlog_vec_sieve(ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order)
{
    int maxtry;
	ulong k, p, p1, pmax, logm1;
	ulong * w;
	dlog_precomp_t pre;
	n_primes_t iter;

    w = flint_malloc( nv * sizeof(ulong));
	for (k = 0; k < nv; k++)
	    w[k] = NOT_FOUND;
	w[1] = 0; 

	/* discrete log on first primes, then sieve */
	pmax = (nv < mod.n) ? nv : mod.n;
	p1 = maxtry = 50; /* FIXME: tune this limit! */
	dlog_precomp_n_init(pre, a, mod.n, na, p1);

	logm1 = (mod.n % 2) ? 0 : dlog_precomp(pre, mod.n - 1);

	n_primes_init(iter);
	while ((p = n_primes_next(iter)) < pmax)
	{ 
        ulong m, wp; 
        if (mod.n % p == 0) /* FIXME: those primes could be known... */
            continue; /* won't be attained another time */
        if (p < p1 || (wp = logp_sieve(w, nv, p, mod.n, logm1, order, maxtry)) != NOT_FOUND)
            wp = nmod_mul(dlog_precomp(pre, p), va, order);
        for (k = p, m = 1; k < nv; k += p, m++)
            w[k] = nmod_add(w[m],  wp, order);
    }
	for (k = 0; k < nv; k++)
        if (v[k] != NOT_FOUND)
            v[k] = nmod_add(v[k], w[k], order);
    n_primes_clear(iter);
    flint_free(w);
}

/* group components up to bound and return cofactor */
#define LOOP 0
#define SIEVE 1
static void
n_factor_group(n_factor_t * fac, ulong bound)
{
    int i, j, k;
    ulong m[FLINT_MAX_FACTORS_IN_LIMB];
    ulong c = 1;
    i = 0;
    for (k = fac->num - 1; k >= 0; k--)
    {
        ulong qe = n_pow(fac->p[k], fac->exp[k]);
        if (qe > bound)
            c *= qe;
        else
        {
            /* try to insert somewhere in m */
            for (j = 0; j < i && (m[j] * qe > bound); j++);
            if (j == i)
                m[i++] = qe;
            else
                m[j] *= qe;
        }
    }
    for (j = 0; j < i; j++)
    {
        fac->p[j] = m[j];
        fac->exp[j] = LOOP;
    }
    if (c > 1)
    {
        fac->p[i] = c;
        fac->exp[i] = SIEVE;
        i++;
    }
    fac->num = i;
}

/* assume v[k] = -1 for bad primes? */
/* loop on small components and keep one subgroup for DLOG + sieve */
  void
dlog_vec_crt(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    n_factor_t fac;
    ulong maxloop;
    int k;

    maxloop = 3 * nv; /* FIXME: tune this */
    n_factor_init(&fac);
    n_factor(&fac, na, 1);
    n_factor_group(&fac, maxloop);
    for (k = 0; k < fac.num; k++)
    {   
        ulong m, M, aM, uM, vaM;
        m = fac.p[k];
        M = na / m;
        aM = nmod_pow_ui(a, M, mod);
        uM = M * n_invmod(M % m, m); /* uM < n */
        vaM = nmod_mul(va, uM % order.n, order);
        if (fac.exp[k] == LOOP)
            dlog_vec_loop(v, nv, aM, vaM, mod, m, order);
        else
            dlog_vec_sieve(v, nv, aM, vaM, mod, m, order);
    }
}
