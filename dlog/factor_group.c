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

#include "dlog.h"

/* group components up to bound and return cofactor */
void
dlog_n_factor_group(n_factor_t * fac, ulong bound)
{
    int i, j, k;
    ulong m[FLINT_MAX_FACTORS_IN_LIMB];
    ulong c = 1;
    i = 0;
    for (k = fac->num - 1; k >= 0; k--)
    {
        ulong qe = n_pow(fac->p[k], fac->exp[k]);
        if (qe > bound)
            c *= qe;
        else
        {
            /* try to insert somewhere in m */
            for (j = 0; j < i && (m[j] * qe > bound); j++);
            if (j == i)
                m[i++] = qe;
            else
                m[j] *= qe;
        }
    }
    for (j = 0; j < i; j++)
    {
        fac->p[j] = m[j];
        fac->exp[j] = G_SMALL;
    }
    if (c > 1)
    {
        fac->p[i] = c;
        fac->exp[i] = G_BIG;
        i++;
    }
    fac->num = i;
}
