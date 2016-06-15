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

#include "acb_dirichlet.h"


/* order of an element knowing the factorization of a multiple */
static ulong
nmod_order_precomp(ulong a, nmod_t mod, ulong expo, n_factor_t fac)
{
    int k;
    ulong pe, ap, order = 1;
    for (k = 0; k < fac.num; k++)
    {
        pe = n_pow(fac.p[k], fac.exp[k]);
        ap = nmod_pow_ui(a, expo / pe, mod);
        while ( ap != 1)
        {
            ap = nmod_pow_ui(ap, fac.p[k], mod);
            order *= fac.p[k];
        }
    }
    return order;
}

static ulong
n_conductor(ulong q, ulong a)
{
    slong k;
    ulong ap, cond;

    nmod_t pe;
    n_factor_t fac;
    n_factor_init(&fac);
    n_factor(&fac, q, 1);

    cond = 1;

    for (k = 0; k < fac.num; k++)
    {
        ulong p, e;
        p = fac.p[k];
        e = fac.exp[k];

        nmod_init(&pe, n_pow(p, e));
        ap = a % pe.n;
        if (ap == 1)
            continue;
        if (p == 2)
        {
            cond = 4;
            if (a % 4 == 3)
                ap = pe.n - ap;
        }
        else
        {
            cond *= p;
            ap = nmod_pow_ui(ap, p - 1, pe);
        }

        while (ap != 1)
        {
            cond *= p;
            ap = nmod_pow_ui(ap, p, pe);
        }

    }

        return cond;
}

int main()
{
    slong iter, bits;
    flint_rand_t state;

    flint_printf("chars....");
    fflush(stdout);
    flint_randinit(state);
    for (bits = 5; bits <= 30; bits += 5)
    { 

    for (iter = 0; iter < 50; iter++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_char_t chi, chi2;
        ulong q, iter2;
        n_factor_t fac;

        q = 2 + n_randint(state, 1 << bits);

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_char_init(chi, G);
        acb_dirichlet_char_init(chi2, G);

        n_factor_init(&fac);
        n_factor(&fac, G->expo, 1);

        acb_dirichlet_group_dlog_precompute(G, 50);

        /* check number char properties */
        for (iter2 = 0; iter2 < 50; iter2++)
        { 
            ulong m, n;
            ulong order, chim1, pairing, cond;

            do
                m = n_randint(state, q);
            while (n_gcd(q, m) > 1);

            acb_dirichlet_char(chi, G, m);

            order = nmod_order_precomp(m, G->mod, G->expo, fac);
            if (order != chi->order)
            {
                flint_printf("FAIL: order\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("order(m) = %wu\n\n", order);
                flint_printf("chi->order = %wu\n\n", chi->order);
                abort();
            }

            cond = n_conductor(G->q, m);
            if (cond != acb_dirichlet_char_conductor(G, chi))
            {
                flint_printf("FAIL: conductor\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("conductor(m) = %wu\n\n", cond);
                flint_printf("chi->conductor = %wu\n\n", acb_dirichlet_char_conductor(G, chi));
                abort();
            }

            chim1 = acb_dirichlet_ui_chi(G, chi, q - 1);
            if (acb_dirichlet_char_parity(chi) != (chim1 != 0))
            {
                flint_printf("FAIL: parity\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("chi(-1) = %wu\n\n", chim1);
                flint_printf("char_parity = %d", acb_dirichlet_char_parity(chi));
                acb_dirichlet_char_print(G, chi);
                abort();
            }

            do
                n = n_randint(state, q);
            while (n_gcd(q, n) > 1);

            acb_dirichlet_char(chi2, G, n);
            pairing = acb_dirichlet_ui_pairing(G, m, n);

            if (pairing != acb_dirichlet_ui_chi(G, chi, n) * (G->expo / chi->order)
                    || pairing != acb_dirichlet_ui_chi(G, chi2, m) * (G->expo / chi2->order))
            {
                flint_printf("FAIL: pairing\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("n = %wu\n\n", n);
                flint_printf("chi(m,n) = %wu\n\n", pairing);
                abort();
            }
        }

        acb_dirichlet_char_clear(chi);
        acb_dirichlet_char_clear(chi2);
        acb_dirichlet_group_clear(G);
    }
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
