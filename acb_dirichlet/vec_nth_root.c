/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* assume xs has size >= len * step */
void
_acb_vec_set_powers_step(acb_ptr xs, slong n, slong len, slong step, slong prec)
{
    slong i, j;
    prec += n_clog(n, 2);

    for (i = 0, j = 0; i < len; i++, j += step)
    {
        if (i == 0)
            acb_one(xs + j);
        else if (i == 1)
            acb_dirichlet_nth_root(xs + j, n, prec);
        else if (i % 2 == 0)
            acb_mul(xs + j, xs + j / 2, xs + j / 2, prec);
        else
            acb_mul(xs + j, xs + j - step, xs + step, prec);
    }
}

/* assume len = p^e and z has size >= len * step */
void
_acb_vec_nth_roots_pe(acb_ptr z, slong p, slong e, slong len, slong step, slong prec)
{
    if (e <= 1)
    {
        _acb_vec_set_powers_step(z, p, len, step, prec);
    }
    else
    {
        slong q, r;
        _acb_vec_nth_roots_pe(z, p, e - 1, len / p, step * p, prec);
        _acb_vec_set_powers_step(z, n_pow(p, e), p, step, prec);

        for (q = p; q < len; q += p)
            for (r = 1; r < p; r++)
                acb_mul(z + (q + r) * step, z + q * step, z + r * step, prec);
    }
}

void
acb_dirichlet_vec_nth_roots(acb_ptr z, slong len, slong prec)
{
   slong i, q;
   n_factor_t fac;

   acb_one(z + 0);

   n_factor_init(&fac);
   n_factor(&fac, len, 0);
   q = 1;

   for (i = 0; i < fac.num; i++)
   {
       slong p, e, pe, mp, mq;
       p = fac.p[i];
       e = fac.exp[i];
       pe = n_pow(p, e);
       mp = len / pe;
       mq = len / q;

       _acb_vec_nth_roots_pe(z, p, e, pe, mp, prec);

       if (i > 0)
       {
           slong k, l;

           for (k = mp; k < len; k += mp)
           {
               for (l = mq; l < len - k; l += mq)
                   acb_mul(z + k + l, z + k, z + l, prec);
               for (; l < len; l += mq)
                   acb_mul(z + k + l - len, z + k, z + l, prec);
           }
       }
       q *= pe;
   }
}
