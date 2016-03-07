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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "arb_mat.h"


int
arb_mat_spd_logdet(arb_t logdet, const arb_mat_t A, slong prec)
{
    slong n;

    if (!arb_mat_is_square(A))
    {
        flint_printf("arb_mat_spd_logdet: a square matrix is required!\n");
        abort();
    }

    if (arb_mat_is_empty(A))
    {
        arb_zero(logdet);
        return 1;
    }

    n = arb_mat_nrows(A);

    if (n == 1)
    {
        if (arb_is_positive(arb_mat_entry(A, 0, 0)))
        {
            arb_log(logdet, arb_mat_entry(A, 0, 0), prec);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        int result;
        arb_mat_t L;
        arb_mat_init(L, n, n);
        result = arb_mat_cho(L, A, prec);
        if (result)
        {
            slong i;
            arb_set(logdet, arb_mat_entry(L, 0, 0));
            for (i = 1; i < n; i++)
            {
                arb_mul(logdet, logdet, arb_mat_entry(L, i, i), prec);
            }
            arb_log(logdet, logdet, prec);
            arb_mul_2exp_si(logdet, logdet, 1);
        }
        arb_mat_clear(L);
        return result;
    }
}
