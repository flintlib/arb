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

#include "double_extras.h"
#include "mag.h"

double
mag_get_d(const mag_t z)
{
    if (mag_is_zero(z))
    {
        return 0.0;
    }
    else if (mag_is_inf(z))
    {
        return D_INF;
    }
    else if (MAG_EXP(z) < -1000 || MAG_EXP(z) > 1000)
    {
        if (fmpz_sgn(MAG_EXPREF(z)) < 0)
            return ldexp(1.0, -1000);
        else
            return D_INF;
    }
    else
    {
        return ldexp(MAG_MAN(z), MAG_EXP(z) - MAG_BITS);
    }
}

