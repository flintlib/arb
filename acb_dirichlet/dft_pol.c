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

#include "acb_poly.h"
#include "acb_dirichlet.h"

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

void
_acb_vec_roots_pe(acb_ptr z, slong p, slong e, slong len, slong step, slong prec)
{
    if (e <= 1)
    {
        _acb_vec_set_powers_step(z, p, len, step, prec);
    }
    else
    {
        slong q, r;
        _acb_vec_roots_pe(z, p, e - 1, len / p, step * p, prec);
        _acb_vec_set_powers_step(z, n_pow(p, e), p, step, prec);

        for (q = p; q < len; q += p)
            for (r = 1; r < p; r++)
                acb_mul(z + (q + r) * step, z + q * step, z + r * step, prec);
    }
}

#if 1
acb_ptr
acb_roots_init(slong len, slong prec)
{
   slong i, q;
   acb_ptr z;
   n_factor_t fac;
   z = _acb_vec_init(len);
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

       _acb_vec_roots_pe(z, p, e, pe, mp, prec);

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

   return z;
}
#else
acb_ptr
acb_roots_init(slong len, slong prec)
{
    acb_t zeta;
    acb_ptr z;
    acb_init(zeta);
    prec += n_clog(len, 2);
    acb_dirichlet_nth_root(zeta, len, prec);
    z = _acb_vec_init(len);
    /* should use factorization */
    _acb_vec_set_powers(z, zeta, len, prec);
    /*
    flint_printf("\nroots [order %ld, prec %ld]\n", len, prec);
    acb_vec_printd(z, len, 30);
    */
    acb_clear(zeta);
    return z;
}
#endif

/* all roots are already computed */
void
_acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, acb_srcptr z, slong len, slong prec)
{
    /* FIXME: huge accuracy loss */
#if 0
    _acb_poly_evaluate_vec_fast(w, v, len, z, len, prec);
#elif 0
    _acb_poly_evaluate_vec_iter(w, v, len, z, len, prec);
#else
    slong i, j;

    for (i = 0; i < len; i++)
    {
        acb_zero(w + i);
        for (j = 0; j < len; j++)
            acb_addmul(w + i, v + j, z + (i * j % len), prec);
    }
#endif
}

void
acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_ptr z;
    z = acb_roots_init(len, prec);
    _acb_dirichlet_dft_pol(w, v, z, len, prec);
    _acb_vec_clear(z, len);
}
