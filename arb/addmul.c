/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_addmul(arb_t c, const arb_t a, const arb_t b)
{
    if (fmpz_is_zero(arb_expref(c)) && fmpz_is_zero(arb_expref(a))
        && fmpz_is_zero(arb_expref(b)))
    {
        /* multiply errors */
        if (fmpz_is_zero(arb_radref(a)))
        {
            if (!fmpz_is_zero(arb_radref(b)))
            {
                /* (c+e) + a * (b+r) = c + a*b + (e+a*r) */
                _fmpz_addmul_abs(arb_radref(c), arb_midref(a), arb_radref(b));
            }
        }
        else
        {
            if (fmpz_is_zero(arb_radref(b)))
            {
                /* (c+e) + (a+r) * b = c + a*b + (e+b*r) */
                _fmpz_addmul_abs(arb_radref(c), arb_midref(b), arb_radref(a));
            }
            else
            {
                /* (a + r) * (b + s) = a*b + a*s + b*r + r*s */
                if (c != a && c != b)
                {
                    /* no aliasing */
                    _fmpz_addmul_abs(arb_radref(c), arb_radref(a), arb_radref(b));
                    _fmpz_addmul_abs(arb_radref(c), arb_midref(a), arb_radref(b));
                    _fmpz_addmul_abs(arb_radref(c), arb_midref(b), arb_radref(a));
                }
                else
                {
                    fmpz_t t;
                    fmpz_init(t);
                    fmpz_mul(t, arb_radref(a), arb_radref(b));
                    _fmpz_addmul_abs(t, arb_midref(a), arb_radref(b));
                    _fmpz_addmul_abs(t, arb_midref(b), arb_radref(a));
                    fmpz_add(arb_radref(c), arb_radref(c), t);
                    fmpz_clear(t);
                }
            }
        }

        /* multiply midpoints */
        fmpz_addmul(arb_midref(c), arb_midref(a), arb_midref(b));

        _arb_normalise(c);
    }
    else
    {
        arb_t t;
        arb_init(t, arb_prec(c));
        arb_mul(t, a, b);
        arb_add(c, c, t);
        arb_clear(t);
    }
}
