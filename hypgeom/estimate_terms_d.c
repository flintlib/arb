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
#include "double_extras.h"
#include "hypgeom.h"

#define LOG2 0.69314718055994530942
#define EXP1 2.7182818284590452354

static __inline__ double d_root(double x, int r)
{
    if (r == 1)
        return x;
    if (r == 2)
        return sqrt(x);
    return pow(x, 1. / r);
}

long
hypgeom_estimate_terms(const mag_t z, int r, slong prec)
{
    double y, t;

    t = mag_get_d(z);

    if (t == 0)
        return 1;

    if (r == 0)
    {
        if (t >= 1)
        {
            flint_printf("z must be smaller than 1\n");
            abort();
        }

        y = (log(1-t) - prec * LOG2) / log(t) + 1;
    }
    else
    {
        y = (prec * LOG2) / (d_root(t, r) * EXP1 * r);
        y = (prec * LOG2) / (r * d_lambertw(y)) + 1;
    }

    if (y >= WORD_MAX / 2)
    {
        flint_printf("error: series will converge too slowly\n");
        abort();
    }

    return y;
}

