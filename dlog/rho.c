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

void
dlog_rho_init(dlog_rho_t t, ulong a, ulong mod, ulong n)
{
  t->a = a;
  t->n = n;
  t->mod = mod;
  t->nisprime = n_is_prime(n);
}

void
dlog_rho_clear(dlog_rho_t t)
{
  return;
}

static ulong
dlog_once(ulong b, ulong a, ulong mod, ulong n)
{
    if (n < 50)
    {
        int k;
        ulong ak = 1;
        for (k=0; k < n; k++)
        {
            if (ak == b)
                return k;
            ak = (ak * a) % mod;
        }
        flint_printf("FAIL[dlog once]: log(%wu,%wu) mod %wu not found (size %wu)\n",
                b, a, mod, n);
        abort();
    } else {
        dlog_rho_t t;
        dlog_rho_init(t, a, mod, n);
        return dlog_rho(t, b);
    }
}

/* solve log knowing equation  e = f * log(b) [n] */
static ulong
dlog_quotient(const dlog_rho_t t, ulong e, ulong f, ulong g, ulong b)
{
    ulong r, n, b_ar, an;
    n = t->n;
    if (g == n)
    {
        flint_printf("FAIL[dlog quotient]: trivial relation e = %wu, f = %wu mod %wu\n",
                e, f, n);
        abort();
    }
    n = n / g;
    e = e / g;
    f = f / g;
    f = n_invmod(f, n);
    r = ( e * f ) % n;
    an = n_powmod(t->a, n, t->mod);
    b_ar = (b * n_invmod(n_powmod(t->a, r, t->mod), t->mod)) % t->mod;
    return r + n * dlog_once(b_ar, an, t->mod, g);
}

#define RWALK 20
ulong
dlog_rho(const dlog_rho_t t, ulong b)
{
    int j, k, l;
    ulong m[RWALK], n[RWALK], ab[RWALK];
    ulong x[2], e[2], f[2], g;
    flint_rand_t state;
    flint_randinit(state);
    do {
        for (k = 0; k < RWALK; k++)
        {
            m[k] = 1 + n_randint(state, t->n - 1);
            n[k] = 1 + n_randint(state, t->n - 1);
            ab[k] = (n_powmod(t->a, m[k], t->mod) * n_powmod(b, n[k], t->mod)) % t->mod;
        }
        /* x[l] = a^e[l] * b^f[l] */
        x[0] = x[1] = 1;
        e[0] = e[1] = 0;
        f[0] = f[1] = 0;
        do {
            for(j = 0; j < 3; j++)
            {
                l = (j > 0);
                k = floor( (double) RWALK * x[l] / t->mod );
                x[l] = (x[l] * ab[k]) % t->mod;
                e[l] = (e[l] + m[k]) % t->n;
                f[l] = (f[l] + n[k]) % t->n;
            }
        } while (x[0] != x[1]);
    } while (e[0] == e[1] && f[0] == f[1]);
    flint_randclear(state);
    /* e = f * log(b) */
    e[0] = (e[0] > e[1]) ?  e[0] - e[1] : e[0] + t->n - e[1];
    f[0] = (f[1] > f[0]) ?  f[1] - f[0] : f[1] + t->n - f[0];
    if (!t->nisprime && (g = n_gcd(f[0], t->n)) > 1)
    {
        return dlog_quotient(t, e[0], f[0], g, b);
    }
    else
    {
        f[0] = n_invmod(f[0], t->n);
        return ( e[0] * f[0] ) % t->n;
    }
}
