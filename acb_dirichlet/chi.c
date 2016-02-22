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

    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_dirichlet.h"

/* todo: modular arithmetic */

static ulong
chi_odd_exponent(const acb_dirichlet_group_t G, ulong m, ulong n)
{
    ulong x, k, pk, gk, logm, logn;

    x = 0;

    for (k = 0; k < G->num; k++)
    {
        pk = n_pow(G->primes[k], G->exponents[k]);
        gk = G->generators[k] % pk;

        logm = n_discrete_log_bsgs(m % pk, gk, pk);
        logn = n_discrete_log_bsgs(n % pk, gk, pk);

        x = (x + G->PHI[k] * logm * logn) % G->phi_q_odd;
    }

    return x;
}

static ulong
chi_even_exponent(const acb_dirichlet_group_t G, ulong m, ulong n)
{
    ulong x;
    ulong q_even = G->q_even;

    if (q_even <= 2)
        return 0;

    x = 0;

    if ((m % 4 == 3) && (n % 4 == 3))
        x = q_even / 8;

    if (q_even > 4)
    {
        ulong g2, logm, logn;

        g2 = 5;

        if (m % 4 == 3)
        {
            m = n_negmod(m, q_even);
        }

        if (n % 4 == 3)
        {
            n = n_negmod(n, q_even);
        }

        logm = n_discrete_log_bsgs(m % q_even, g2, q_even);
        logn = n_discrete_log_bsgs(n % q_even, g2, q_even);
        x += logm * logn;
    }

    return x % (q_even / 4);
}

void
acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec)
{
    fmpq_t t, u;

    ulong odd_part, even_part;
    ulong q_even = G->q_even;
    ulong q_odd = G->q_odd;

    if ((q_even > 1 && (n % 2 == 0)) || (q_odd > 1 && (n_gcd(q_odd, n) != 1)))
    {
        acb_zero(res);
        return;
    }

    if (q_even > 2)
    {
        if (q_even == 4)
        {
            if (m % 4 == 3 && n % 4 == 3)
                even_part = q_even / 2;  /* -1 */
            else
                even_part = 0;           /* 1 */
        }
        else
        {
            even_part = 4 * chi_even_exponent(G, m % q_even, n % q_even);
        }
    }
    else
    {
        even_part = 0;
    }

    if (q_odd > 1)
        odd_part = chi_odd_exponent(G, m % q_odd, n % q_odd);
    else
        odd_part = 0;

    fmpq_init(t);
    fmpq_init(u);

    fmpq_set_si(t, 2 * even_part, q_even);
    fmpq_set_si(u, 2 * odd_part, G->phi_q_odd);
    fmpq_add(t, t, u);

    arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), t, prec);

    fmpq_clear(t);
    fmpq_clear(u);
}

