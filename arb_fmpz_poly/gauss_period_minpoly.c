/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "acb_dirichlet.h"

void
arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n)
{
    ulong k, d, e, g, gk, qinv;
    ulong * es;
    slong prec, initial_prec;
    int done, real;
    int * lower_plane;

    if (n == 0 || !n_is_prime(q) || ((q - 1) % n) != 0 ||
            n_gcd_full(n, (q - 1) / n) != 1 || q >= ULONG_MAX / 2)
    {
        fmpz_poly_zero(res);
        return;
    }

    d = (q - 1) / n;

    /* this is much faster */
    if (d == 1)
    {
        fmpz_poly_cyclotomic(res, q);
        return;
    }

    g = n_primitive_root_prime(q);
    qinv = n_preinvert_limb(2 * q);

    es = flint_malloc(sizeof(ulong) * d);
    lower_plane = flint_calloc(n, sizeof(int));

    for (e = 0; e < d; e++)
        es[e] = n_powmod2(g, n * e, 2 * q);

    /* either all roots are real, or all roots are complex */
    real = (n % 2) == 1;

    /* first estimate precision crudely based on d and n */
    initial_prec = n * log(2 * d) * 1.4426950408889 * 1.03 + 20;
    initial_prec = FLINT_MAX(initial_prec, 48);

    /* if high, start lower to get a good estimate */
    if (initial_prec > 200)
        initial_prec = 48;

    for (prec = initial_prec, done = 0; !done; )
    {
        acb_dirichlet_roots_t zeta;
        arb_poly_t pz;
        arb_ptr roots;
        acb_ptr croots;
        acb_t t, u;
        slong root_index;

        acb_dirichlet_roots_init(zeta, q, n * d, prec);
        roots = _arb_vec_init(n);
        croots = (acb_ptr) roots;

        acb_init(t);
        acb_init(u);
        arb_poly_init(pz);

        root_index = 0;

        for (k = 0; k < n; k++)
        {
            if (lower_plane[k])
                continue;

            gk = n_powmod2(g, k, 2 * q);
            acb_zero(u);

            if (real)
            {
                for (e = 0; e < d / 2; e++)
                {
                    acb_dirichlet_root(t, zeta, n_mulmod2_preinv(gk, es[e], 2 * q, qinv), prec);
                    acb_mul_2exp_si(t, t, 1);  /* compute conjugates */
                    acb_add(u, u, t, prec);
                }

                arb_set(roots + k, acb_realref(u));
            }
            else
            {
                for (e = 0; e < d; e++)
                {
                    acb_dirichlet_root(t, zeta, n_mulmod2_preinv(gk, es[e], 2 * q, qinv), prec);
                    acb_add(u, u, t, prec);
                }

                if (arb_is_negative(acb_imagref(u)))
                {
                    lower_plane[k] = 1;
                }
                else if (arb_contains_zero(acb_imagref(u)))
                {
                    /* todo: could increase precision */
                    flint_printf("fail! imaginary part should be nonzero\n");
                    flint_abort();
                }
                else
                {
                    acb_set(croots + root_index, u);
                    root_index++;
                }
            }
        }

        if (real)
            arb_poly_product_roots(pz, roots, n, prec);
        else
            arb_poly_product_roots_complex(pz, NULL, 0, croots, root_index, prec);

        done = arb_poly_get_unique_fmpz_poly(res, pz);

        if (!done && prec == initial_prec)
        {
            mag_t m, mmax;
            mag_init(m);
            mag_init(mmax);

            for (k = 0; k < n; k++)
            {
                arb_get_mag(m, pz->coeffs + k);
                mag_max(mmax, mmax, m);
            }

            prec = mag_get_d_log2_approx(mmax) * 1.03 + 20;

            if (prec < 2 * initial_prec)
                prec = 2 * initial_prec;

            mag_clear(m);
            mag_clear(mmax);
        }
        else if (!done)
        {
            prec *= 2;
        }

        acb_dirichlet_roots_clear(zeta);
        _arb_vec_clear(roots, n);
        acb_clear(t);
        acb_clear(u);
        arb_poly_clear(pz);
    }

    flint_free(es);
    flint_free(lower_plane);
}

