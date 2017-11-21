/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"
#include "acb_calc.h"

/*
  Gauss-Legendre quadrature nodes are cached to speed up multiple integrations
  and adaptive subdivision. The steps of 2^(n/2) used here give slightly better
  performance than steps of 2^n (we use at most 1.4x more points than
  needed and not 2x more points) but may require more precomputation.
*/

#define GL_STEPS 38

const slong gl_steps[GL_STEPS] = {1, 2, 4, 6, 8, 12, 16, 22, 32, 46, 64,
    90, 128, 182, 256, 362, 512, 724, 1024, 1448, 2048, 2896, 4096,
    5792, 8192, 11586, 16384, 23170, 32768, 46340, 65536, 92682,
    131072, 185364, 262144, 370728, 524288, 741456};

typedef struct
{
    slong gl_prec[GL_STEPS];
    arb_ptr gl_nodes[GL_STEPS];
    arb_ptr gl_weights[GL_STEPS];
}
gl_cache_struct;

TLS_PREFIX gl_cache_struct * gl_cache = NULL;

void gl_cleanup()
{
    slong i;

    if (gl_cache == NULL)
        return;

    for (i = 0; i < GL_STEPS; i++)
    {
        if (gl_cache->gl_prec[i] != 0)
        {
            _arb_vec_clear(gl_cache->gl_nodes[i], (gl_steps[i] + 1) / 2);
            _arb_vec_clear(gl_cache->gl_weights[i], (gl_steps[i] + 1) / 2);
        }
    }

    flint_free(gl_cache);
    gl_cache = NULL;
}

void gl_init()
{
    gl_cache = flint_calloc(1, sizeof(gl_cache_struct));
    flint_register_cleanup_function(gl_cleanup);
}

/* Compute GL node and weight of index k for n = gl_steps[i]. Cached. */
void
acb_calc_gl_node(arb_t x, arb_t w, slong i, slong k, slong prec)
{
    slong n, kk, jj, wp;

    if (i < 0 || i >= GL_STEPS || prec < 2)
        flint_abort();

    if (gl_cache == NULL)
        gl_init();

    n = gl_steps[i];

    if (k < 0 || k >= n)
        flint_abort();

    if (2 * k < n)
        kk = k;
    else
        kk = n - 1 - k;

    if (gl_cache->gl_prec[i] < prec)
    {
        if (gl_cache->gl_prec[i] == 0)
        {
            gl_cache->gl_nodes[i] = _arb_vec_init((n + 1) / 2);
            gl_cache->gl_weights[i] = _arb_vec_init((n + 1) / 2);
        }

        wp = FLINT_MAX(prec, gl_cache->gl_prec[i] * 2 + 30);

        for (jj = 0; 2 * jj < n; jj++)
        {
            arb_hypgeom_legendre_p_ui_root(gl_cache->gl_nodes[i] + jj,
                gl_cache->gl_weights[i] + jj, n, jj, wp);
        }

        gl_cache->gl_prec[i] = wp;
    }

    if (2 * k < n)
        arb_set_round(x, gl_cache->gl_nodes[i] + kk, prec);
    else
        arb_neg_round(x, gl_cache->gl_nodes[i] + kk, prec);

    arb_set_round(w, gl_cache->gl_weights[i] + kk, prec);
}

int
acb_calc_integrate_gl_auto_deg(acb_t res, slong * eval_count,
    acb_calc_func_t f, void * param,
    const acb_t a, const acb_t b, const mag_t tol,
    slong deg_limit, int verbose, slong prec)
{
    acb_t mid, delta, wide;
    mag_t tmpm;
    slong status;
    acb_t s, v;
    mag_t M, X, Y, rho, err, t, best_rho;
    slong k, Xexp;
    slong i, n, best_n;

    status = ARB_CALC_NO_CONVERGENCE;

    if (deg_limit <= 0)
    {
        acb_indeterminate(res);
        eval_count[0] = 0;
        return status;
    }

    acb_init(mid);
    acb_init(delta);
    acb_init(wide);
    mag_init(tmpm);

    /* delta = (b-a)/2 */
    acb_sub(delta, b, a, prec);
    acb_mul_2exp_si(delta, delta, -1);

    /* mid = (a+b)/2 */
    acb_add(mid, a, b, prec);
    acb_mul_2exp_si(mid, mid, -1);

    acb_init(s);
    acb_init(v);
    mag_init(M);
    mag_init(X);
    mag_init(Y);
    mag_init(rho);
    mag_init(t);
    mag_init(err);
    mag_init(best_rho);

    best_n = -1;
    eval_count[0] = 0;

    mag_inf(err);

    for (Xexp = 0; Xexp < prec /* && Xexp == 0 */; Xexp += FLINT_MAX(1, Xexp))
    {
        mag_one(X);
        mag_mul_2exp_si(X, X, Xexp + 1);

        /* rho = X + sqrt(X^2 - 1)  (lower bound) */
        mag_mul_lower(rho, X, X);
        mag_one(t);
        mag_sub_lower(rho, rho, t);
        mag_sqrt_lower(rho, rho);
        mag_add_lower(rho, rho, X);

        /* Y = sqrt(X^2 - 1)  (upper bound) */
        mag_mul(Y, X, X);
        mag_one(t);
        mag_sub(Y, Y, t);
        mag_sqrt(Y, Y);

        acb_zero(wide);
        mag_set(arb_radref(acb_realref(wide)), X);
        mag_set(arb_radref(acb_imagref(wide)), Y);

        /* transform to [a,b] */
        acb_mul(wide, wide, delta, prec);
        acb_add(wide, wide, mid, prec);

        f(v, wide, param, 1, prec);
        eval_count[0]++;

        /* no chance */
        if (!acb_is_finite(v))
            break;

        /* M = (b-a)/2  |f| */
        acb_get_mag(M, v);
        acb_get_mag(tmpm, delta);
        mag_mul(M, M, tmpm);

        /* Search for the smallest n that gives err < tol (if possible) */
        for (i = 0; i < GL_STEPS && gl_steps[i] <= deg_limit; i++)
        {
            n = gl_steps[i];

            /* (64/15) M / ((rho-1) rho^(2n-1)) */
            mag_pow_ui_lower(t, rho, 2 * n - 1);
            mag_one(tmpm);
            mag_sub_lower(tmpm, rho, tmpm);
            mag_mul_lower(t, t, tmpm);
            mag_mul_ui_lower(t, t, 15);
            mag_div(t, M, t);
            mag_mul_2exp_si(t, t, 6);

            if (mag_cmp(t, tol) < 0)
            {
                status = ARB_CALC_SUCCESS;

                /* The best so far. */
                if (best_n == -1 || n < best_n)
                {
                    mag_set(err, t);
                    mag_set(best_rho, rho);
                    best_n = n;
                }

                /* Best possible n. */
                if (n == 1)
                    break;
            }
        }
    }

    /* Evaluate best found Gauss-Legendre quadrature rule. */
    if (status == ARB_CALC_SUCCESS)
    {
        arb_t x, w;
        arb_init(x);
        arb_init(w);

        if (verbose)
        {
            acb_get_mag(tmpm, delta);
            flint_printf("  {GL deg %ld on [", best_n);
            acb_printn(a, 10, ARB_STR_NO_RADIUS); flint_printf(", ");
            acb_printn(b, 10, ARB_STR_NO_RADIUS);
            flint_printf("], delta "); mag_printd(tmpm, 5);
            flint_printf(", rho "); mag_printd(best_rho, 5);
            flint_printf(", tol "); mag_printd(tol, 3);
            flint_printf("}\n");
        }

        if (best_n == -1)
            flint_abort();

        for (i = 0; i < GL_STEPS; i++)
            if (gl_steps[i] == best_n)
                break;

        acb_zero(s);

        for (k = 0; k < best_n; k++)
        {
            acb_calc_gl_node(x, w, i, k, prec);
            acb_mul_arb(wide, delta, x, prec);
            acb_add(wide, wide, mid, prec);
            f(v, wide, param, 0, prec);
            acb_addmul_arb(s, v, w, prec);
        }

        eval_count[0] += best_n;

        acb_mul(res, s, delta, prec);
        acb_add_error_mag(res, err);

        arb_clear(x);
        arb_clear(w);
    }
    else
    {
        acb_indeterminate(res);
    }

    acb_clear(s);
    acb_clear(v);
    mag_clear(M);
    mag_clear(X);
    mag_clear(Y);
    mag_clear(rho);
    mag_clear(t);
    mag_clear(err);
    mag_clear(best_rho);

    acb_clear(mid);
    acb_clear(delta);
    acb_clear(wide);
    mag_clear(tmpm);

    return status;
}

