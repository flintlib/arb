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
dlog_power_init(dlog_power_t t, ulong a, ulong mod, ulong p, ulong e, ulong num)
{
    int k;
    nmod_init(&t->mod, mod);
    t->p = p;
    t->e = e;

    t->apk = flint_malloc(e * sizeof(ulong));
    t->pre = flint_malloc(sizeof(dlog_precomp_struct));

    t->apk[0] = nmod_inv(a, t->mod);
    for (k = 1; k < e; k++)
        t->apk[k] = nmod_pow_ui(t->apk[k-1], p, t->mod);

    dlog_precomp_p_init(t->pre, nmod_inv(t->apk[e-1], t->mod), mod, p, num);
}

void
dlog_power_clear(dlog_power_t t)
{
    flint_free(t->apk);
    dlog_precomp_clear(t->pre);
    flint_free(t->pre);
}

ulong
dlog_power(const dlog_power_t t, ulong b)
{
    int k;
    ulong x, pk[30]; /* 3^30*2+1, 2^30*3+1 are primes */

    pk[0] = 1;
    for (k = 1; k < t->e; k++)
       pk[k] = pk[k-1] * t->p;
    
    x = 0;
    for(k = 0; k < t->e; k++)
    {
      ulong bk, xk;
      bk = nmod_pow_ui(b, pk[t->e-1-k], t->mod);
      xk = dlog_precomp(t->pre, bk);
      b = nmod_mul(b, nmod_pow_ui(t->apk[k], xk, t->mod), t->mod);
      x += xk * pk[k]; /* cannot overflow */
    }
    return x;
}
