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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include "arith.h"
#include "fmprb.h"

void
fmprb_zeta_ui(fmprb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmprb_set_si(x, -1);
        fmprb_mul_2exp_si(x, x, -1);
    }
    else if (n == 1)
    {
        printf("exception: zeta_ui(1)\n");
        abort();
    }
    /* fast detection of asymptotic case */
    else if (n > 0.7 * prec)
    {
        fmprb_zeta_ui_asymp(x, n, prec);
    }
    else if (n == 3)
    {
        fmprb_const_zeta3_bsplit(x, prec);
    }
    /* small even n */
    else if ((n % 2 == 0) && (n < 40 + 0.11*prec))
    {
        fmprb_zeta_ui_bernoulli(x, n, prec);
    }
    /* small odd n, extremely high precision */
    else if (n < prec * 0.0006)
    {
        fmprb_zeta_ui_bsplit(x, n, prec);
    }
    /* large n */
    else if (prec > 20 && n > 6 && n > 0.4 * pow(prec, 0.8))
    {
        fmprb_zeta_ui_euler_product(x, n, prec);
    }
    /* fallback */
    else
    {
        fmprb_zeta_ui_vec_borwein(x, n, 1, 0, prec);
    }
}

void
fmprb_zeta_ui_vec_even(fmprb_struct * x, ulong start, long num, long prec)
{
    long i;

    for (i = 0; i < num; i++)
        fmprb_zeta_ui(x + i, start + 2 * i, prec);
}

void
fmprb_zeta_ui_vec_odd(fmprb_struct * x, ulong start, long num, long prec)
{
    long i, num_borwein;
    ulong cutoff;

    cutoff = 40 + 0.3 * prec;

    if (cutoff > start)
    {
        num_borwein = 1 + (cutoff - start) / 2;
        num_borwein = FLINT_MIN(num_borwein, num);
    }
    else
        num_borwein = 0;

    fmprb_zeta_ui_vec_borwein(x, start, num_borwein, 2, prec);
    for (i = num_borwein; i < num; i++)
        fmprb_zeta_ui(x + i, start + 2 * i, prec);
}

void
fmprb_zeta_ui_vec(fmprb_struct * x, ulong start, long num, long prec)
{
    long i, num_odd, num_even, start_odd;
    fmprb_struct * tmp;

    num_odd = num / 2 + (start & num & 1);
    num_even = num - num_odd;

    start_odd = start % 2;

    fmprb_zeta_ui_vec_even(x, start + start_odd, num_even, prec);
    fmprb_zeta_ui_vec_odd(x + num_even, start + !start_odd, num_odd, prec);

    /* interleave */
    tmp = flint_malloc(sizeof(fmprb_struct) * num);
    for (i = 0; i < num_even; i++) tmp[i] = x[i];
    for (i = 0; i < num_odd; i++) tmp[num_even + i] = x[num_even + i];
    for (i = 0; i < num_even; i++) x[start_odd + 2  * i] = tmp[i];
    for (i = 0; i < num_odd; i++) x[!start_odd + 2  * i] = tmp[num_even + i];
    flint_free(tmp);
}
