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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"

/* When splitting [a,b) into [a,m), [m,b), we need the power x^(m-a).
   This function computes all the exponents (m-a) that can appear when
   doing binary splitting with the top-level interval [0,n),
   assuming that we always choose m = a + floor((b-a)/2), and that
   the case b-a = 2 is inlined. */
long
_arb_compute_bs_exponents(long * tab, long n)
{
    long a, b, aa, ab, ba, bb, length;

    if (n == 1)
    {
        tab[0] = 1;
        return 1;
    }

    if (n == 2 || n == 3 || n == 4)
    {
        tab[0] = 1;
        tab[1] = 2;
        return 2;
    }

    if (n == 6)
    {
        tab[0] = 1;
        tab[1] = 2;
        tab[2] = 3;
        return 3;
    }

    /* first binary splitting call */
    a = n >> 1;
    b = n - (n >> 1);
    tab[0] = a;
    length = 1;

    for (;;)
    {
        /* split a -> aa, ab  and b -> ba, bb */
        aa = a >> 1;
        ab = a - aa;
        ba = b >> 1;
        bb = b - ba;

        tab[length] = ba;
        length++;

        /* at length 3, we split into 2, 1 (and maybe also 2, 2) */
        if (ba == 3)
        {
            tab[length] = 2;
            tab[length + 1] = 1;
            length += 2;
            break;
        }

        /* stop if we have reached 1, or if at 2 and the
           length is a power of 2 (in which case we never reach 1) */
        if (ba == 1 || (ba == 2 && (n & (n-1)) == 0))
            break;

        /* if left and right lengths are different, also add the
           right length */
        if (aa != ba && aa != 1)
        {
            tab[length] = aa;
            length++;
        }

        a = aa;
        b = bb;
    }

    /* Prevent compiler warning (ab is used to make the algorithm clear). */
    ab = ab;

    /* we always include x^1 in the table, even if the binary splitting
       terminates at step length 2 */
    if (tab[length-1] != 1)
    {
        tab[length] = 1;
        length++;
    }

    /* reverse table */
    for (a = 0; a < length / 2; a++)
    {
        b = tab[a];
        tab[a] = tab[length - a - 1];
        tab[length - a - 1] = b;
    }

    return length;
}

/* just do a linear search */
long
_arb_get_exp_pos(const long * tab, long step)
{
    long i;

    for (i = 0; ; i++)
    {
        if (tab[i] == step)
            return i;

        if (tab[i] == 0)
        {
            printf("ERROR: exponent %ld not in table!\n", step);
            abort();
        }
    }
}

static void
bsplit(fmpz_t T, fmpz_t Q, mp_bitcnt_t * Qexp,
    const long * xexp,
    const fmpz * xpow, mp_bitcnt_t r, long a, long b)
{
    int cc;

    if (b - a == 1)
    {
        count_trailing_zeros(cc, (a + 1));
        fmpz_set_ui(Q, (a + 1) >> cc);
        *Qexp = r + cc;

        fmpz_set(T, xpow);
    }
    else if (b - a == 2)
    {
        fmpz_mul_ui(T, xpow, a + 2);
        fmpz_mul_2exp(T, T, r);
        fmpz_add(T, T, xpow + 1);

        count_trailing_zeros(cc, (a + 2));
        fmpz_set_ui(Q, (a + 2) >> cc);
        *Qexp = r + cc;

        count_trailing_zeros(cc, (a + 1));
        fmpz_mul_ui(Q, Q, (a + 1) >> cc);
        *Qexp += r + cc;
    }
    else
    {
        long step, m, i;
        mp_bitcnt_t Q2exp[1];
        fmpz_t Q2, T2;

        step = (b - a) / 2;
        m = a + step;

        fmpz_init(Q2);
        fmpz_init(T2);

        bsplit(T,  Q,  Qexp,  xexp, xpow, r, a, m);
        bsplit(T2, Q2, Q2exp, xexp, xpow, r, m, b);

        fmpz_mul(T, T, Q2);
        fmpz_mul_2exp(T, T, *Q2exp);

        /* find x^step in table */
        i = _arb_get_exp_pos(xexp, step);
        fmpz_addmul(T, xpow + i, T2);  
        fmpz_clear(T2);

        fmpz_mul(Q, Q, Q2);
        *Qexp = *Qexp + *Q2exp;
        fmpz_clear(Q2);
    }
}

void
_arb_exp_sum_bs_powtab(fmpz_t T, fmpz_t Q, mp_bitcnt_t * Qexp,
    const fmpz_t x, mp_bitcnt_t r, long N)
{
    long * xexp;
    long length, i;
    fmpz * xpow;

    /* compute the powers of x that will appear (at least x^1) */
    xexp = flint_calloc(2 * FLINT_BITS, sizeof(long));
    length = _arb_compute_bs_exponents(xexp, N);

    xpow = _fmpz_vec_init(length);
    xpow[0] = *x;   /* create shallow copy of x */

    /* build x^i table */
    for (i = 1; i < length; i++)
    {
        if (xexp[i] == 2 * xexp[i-1])
        {
            fmpz_mul(xpow + i, xpow + i - 1, xpow + i - 1);
        }
        else if (xexp[i] == 2 * xexp[i-2]) /* prefer squaring if possible */
        {
            fmpz_mul(xpow + i, xpow + i - 2, xpow + i - 2);
        }
        else if (xexp[i] == 2 * xexp[i-1] + 1)
        {
            fmpz_mul(xpow + i, xpow + i - 1, xpow + i - 1);
            fmpz_mul(xpow + i, xpow + i, xpow);
        }
        else if (xexp[i] == 2 * xexp[i-2] + 1)
        {
            fmpz_mul(xpow + i, xpow + i - 2, xpow + i - 2);
            fmpz_mul(xpow + i, xpow + i, xpow);
        }
        else
        {
            printf("power table has the wrong structure!\n");
            abort();
        }
    }

    bsplit(T, Q, Qexp, xexp, xpow, r, 0, N);

    fmpz_init(xpow + 0);  /* don't free the shallow copy of x */
    _fmpz_vec_clear(xpow, length);
    flint_free(xexp);
}

