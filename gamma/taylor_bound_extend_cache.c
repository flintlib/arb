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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "gamma.h"

__thread long * gamma_taylor_bound_mag_cache = NULL;

__thread fmpr_struct * gamma_taylor_bound_ratio_cache = NULL;

__thread long gamma_taylor_bound_cache_num = 0;

void
gamma_taylor_bound_extend_cache(long n)
{
    if (n >= gamma_taylor_bound_cache_num)
    {
        long num, i;

        num = FLINT_MAX(n + 1, 2 * gamma_taylor_bound_cache_num);
        num = FLINT_MAX(num, 8);

        gamma_taylor_bound_mag_cache = flint_realloc(
            gamma_taylor_bound_mag_cache, num * sizeof(long));
        gamma_taylor_bound_ratio_cache = flint_realloc(
            gamma_taylor_bound_ratio_cache, num * sizeof(fmpr_struct));

        for (i = gamma_taylor_bound_cache_num; i < num; i++)
        {
            fmpr_init(gamma_taylor_bound_ratio_cache + i);

            if (i > 0)
            {
                gamma_taylor_bound_mag_cache[i] = gamma_taylor_bound_mag(i);
                gamma_taylor_bound_ratio(gamma_taylor_bound_ratio_cache + i, i);
            }
        }

        gamma_taylor_bound_cache_num = num;
    }
}

