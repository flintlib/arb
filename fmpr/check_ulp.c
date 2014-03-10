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

#include "fmpr.h"

int
fmpr_check_ulp(const fmpr_t result, long r, long prec)
{
    if (r == FMPR_RESULT_EXACT)
    {
        return 1;
    }
    else
    {

        fmpr_t err;
        fmpr_t ulp;
        int ok;

        fmpr_init(err);
        fmpr_init(ulp);

        fmpr_ulp(ulp, result, prec);
        fmpr_set_error_result(err, result, r);

        ok = fmpr_equal(err, ulp);

        fmpr_clear(err);
        fmpr_clear(ulp);

        return ok;
    }
}

