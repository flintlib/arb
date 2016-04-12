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
#include "acb_poly.h"

/* x = Pi / q * t^2 */
static void
acb_dirichlet_arb_theta_argt(arb_t x, ulong q, const arb_t t, slong prec)
{
    arb_const_pi(x, prec);
    arb_div_ui(x, x, q, prec);
    arb_mul(x, x, t, prec);
    arb_mul(x, x, t, prec);
    arb_neg(x, x);
    arb_exp(x, x, prec);
}

/* small order, multiply by chi at the end */
void
acb_dirichlet_arb_theta_smallorder(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_powers_t z, slong len, slong prec)
{
    slong k;
    ulong order = z->order;
    arb_t xk2, kxk2, dx, x2;
    acb_ptr t;
    arb_init(xk2);
    arb_init(dx);
    arb_init(x2);
    arb_init(kxk2);

    arb_set(dx, x);
    arb_set(xk2, x);
    arb_mul(x2, x, x, prec);

    t = _acb_vec_init(order);
    _acb_vec_zero(t, order);

    arb_set(acb_realref(t + 0), xk2);

    /* TODO: reduce precision at each step */
    for (k = 2; k < len; k++)
    {
        arb_mul(dx, dx, x2, prec);
        arb_mul(xk2, xk2, dx, prec);
        if (a[k] != ACB_DIRICHLET_CHI_NULL)
        {
           if (parity)
           {
               arb_mul_si(kxk2, xk2, k, prec);
               arb_add(acb_realref(t + a[k]), acb_realref(t + a[k]), kxk2, prec);
           }
           else
           {
               arb_add(acb_realref(t + a[k]), acb_realref(t + a[k]), xk2, prec);
           }
        }
    }

    /* now Hörner */
    _acb_poly_evaluate(res, t, order, z->z + 1, prec);

    _acb_vec_clear(t, order);
    arb_clear(xk2);
    arb_clear(x2);
    arb_clear(dx);
}

void
acb_dirichlet_arb_theta_naive(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_powers_t z, slong len, slong prec)
{
    slong k;
    arb_t xk2, dx, x2;
    acb_t zk;
    arb_init(xk2);
    arb_init(dx);
    arb_init(x2);
    acb_init(zk);

    arb_set(dx, x);
    arb_set(xk2, dx);
    arb_mul(x2, dx, dx, prec);

    acb_set_arb(res, xk2);

    /* TODO: reduce prec */
    for (k = 2; k < len; k++)
    {
        arb_mul(dx, dx, x2, prec);
        arb_mul(xk2, xk2, dx, prec);
        if (a[k] != ACB_DIRICHLET_CHI_NULL)
        {
            acb_dirichlet_power(zk, z, a[k], prec); 
            if (parity)
                acb_mul_si(zk, zk, k, prec);
            acb_addmul_arb(res, zk, xk2, prec);
        }
    }

    arb_clear(xk2);
    arb_clear(x2);
    arb_clear(dx);
    acb_clear(zk);
}

void
acb_dirichlet_chi_theta(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t t, slong prec)
{
    slong len;
    ulong * a;
    arb_t x;
    acb_dirichlet_powers_t z;

    len = acb_dirichlet_theta_length(G->q, t, prec);

    a = flint_malloc(len * sizeof(ulong));
    acb_dirichlet_ui_chi_vec(a, G, chi, len);
    acb_dirichlet_powers_init(z, chi->order, len, prec);

    arb_init(x);
    acb_dirichlet_arb_theta_argt(x, G->q, t, prec);
    acb_dirichlet_arb_theta_naive(res, x, chi->parity, a, z, len, prec);

    arb_clear(x);
    flint_free(a);
    acb_dirichlet_powers_clear(z);
}
