/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "acb_dirichlet.h"

static int usage(char *argv[])
{
    printf("usage: %s [--progress <step>] [--all] [--odd] [--prec <bits>] [--onlymod] {qmin qmax | --value q a}\n", argv[0]);
    return 1;
}

static void
value(ulong q, ulong a, slong prec, slong digits)
{
    acb_dirichlet_group_t G;
    acb_dirichlet_char_t chi;
    acb_t res;
    arb_t one;

    acb_dirichlet_group_init(G, q);
    acb_dirichlet_char_init(chi, G);
    acb_dirichlet_char(chi, G, a);
    acb_init(res);
    arb_init(one);

    arb_one(one);
    acb_dirichlet_theta_arb(res, G, chi, one, prec);

    acb_printd(res, digits);
    flint_printf("\n");

    arb_clear(one);
    acb_clear(res);
    acb_dirichlet_group_clear(G);
    acb_dirichlet_char_clear(chi);
}

static void
check_q(ulong q, int odd, slong prec, slong digits, int onlymod)
{
    slong s;
    ulong k, len;
    acb_dirichlet_group_t G;
    acb_dirichlet_conrey_t x;
    acb_dirichlet_group_init(G, q);
    acb_dirichlet_conrey_init(x, G);
    acb_ptr theta, z;
    arb_t z1;
    arb_ptr t;

    theta = _acb_vec_init(G->phi_q);
    z =  _acb_vec_init(G->phi_q);

    arb_init(z1);
    arb_const_pi(z1, prec);
    arb_div_ui(z1, z1, q, prec);
    arb_neg(z1, z1);
    arb_exp(z1, z1, prec);

    /* len = acb_dirichlet_theta_length_d(q, 1., prec); */
    t = _arb_vec_init(q);
    acb_dirichlet_arb_quadratic_powers(t, q, z1, prec);

    for (s = 0; s <= odd; s++)
    {
        k = 0;
        acb_dirichlet_conrey_one(x, G);
        do {
            acb_set_arb(z + k, t + x->n);
            k++;
        } while (acb_dirichlet_conrey_next(x, G) >= 0);

        acb_dirichlet_dft_conrey(theta, z, G, prec);

        /* check zeros */
        acb_dirichlet_conrey_one(x, G);
        for (k = 0; k < G->phi_q; k++)
        {
            if (acb_contains_zero(theta + k)
                    && acb_dirichlet_conrey_conductor(G, x) == q
                    && acb_dirichlet_conrey_parity(G, x) == s)
            {
                if (onlymod)
                {
                    flint_printf("%wu,%wu\n",q,x->n);
                }
                else
                {
                    flint_printf("\ntheta null q = %wu, n = %wu\n",q, x->n);
                    acb_printd(theta + k, digits);
                    flint_printf("\n");
                }
            }
            acb_dirichlet_conrey_next(x, G);
        }

        if (odd)
        {
            /* change t for odd characters */
            for (k = 0; k < q; k++)
                arb_mul_ui(t + k, t + k, k, prec);
        }

    }

    _arb_vec_clear(t, q);
    _acb_vec_clear(theta, G->phi_q);
    _acb_vec_clear(z, G->phi_q);
    arb_clear(z1);
    acb_dirichlet_conrey_clear(x);
    acb_dirichlet_group_clear(G);
}


int main(int argc, char *argv[])
{
    int i, all = 0, odd = 0, onlymod = 0, eval = 0;
    slong prec = 50, digits = 10;
    slong step = 0;
    ulong q, qmin, qmax;
    n_primes_t iter;

    if (argc < 3)
        return usage(argv);

    for (i = 1; i < argc - 2; i++)
    {
        if (!strcmp(argv[i],"--eval"))
            eval = 1;
        if (!strcmp(argv[i],"--all"))
            all = 1;
        if (!strcmp(argv[i],"--onlymod"))
            onlymod = 1;
        else if (!strcmp(argv[i],"--odd"))
            odd = 1;
        else if (!strcmp(argv[i],"--progress"))
        {
            i++;
            step = atol(argv[i]);
        }
        else if (!strcmp(argv[i],"--prec"))
        {
            i++;
            prec = atol(argv[i]);
            digits = floor(prec * 0.3);
        }
    }

    if (argc < i + 2)
        return usage(argv);

    qmin = atol(argv[i]);
    qmax = atol(argv[i+1]);

    if (eval)
    {
        value(qmin, qmax, prec, digits);
    }
    else
    {
        if (all)
        {
            for (q = qmin; q <= qmax; q++)
            {
                if (q % 4 == 2)
                    continue;
                check_q(q, odd, prec, digits, onlymod);
                if (step && q % step == 0)
                    flint_printf("[%wu]",q);
            }
        }
        else
        {
            ulong p;
            slong it = 0;
            /* look for vanishing theta values for prime power moduli */
            n_primes_init(iter);

            while ((p = n_primes_next(iter)) < qmax)
            {

                for (q = p; q < qmin; q*= p);
                for (; q < qmax; q *= p)
                    check_q(q, odd, prec, digits, onlymod);

                if (step && (it++ % step == 0))
                    flint_printf("[%wu]",p);
            }

            n_primes_clear(iter);
        }
    }

    flint_cleanup();
    return EXIT_SUCCESS;
}
