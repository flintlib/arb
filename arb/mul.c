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
arb_mul(arb_t c, const arb_t a, const arb_t b)
{
    /* multiply errors */
    if (fmpz_is_zero(arb_radref(a)))
    {
        if (fmpz_is_zero(arb_radref(b)))
        {
            fmpz_zero(arb_radref(c));
        }
        else
        {
            /* a * (b + r) = a*b + a*r */
            _fmpz_mul_abs(arb_radref(c), arb_midref(a), arb_radref(b));
        }
    }
    else
    {
        if (fmpz_is_zero(arb_radref(b)))
        {
            /* (a + r) * b = a*b + b*r */
            _fmpz_mul_abs(arb_radref(c), arb_midref(b), arb_radref(a));
        }
        else
        {
            /* (a + r) * (b + s) = a*b + a*s + b*r + r*s */
            if (c != a && c != b)
            {
                /* no aliasing */
                fmpz_mul(arb_radref(c), arb_radref(a), arb_radref(b));
                _fmpz_addmul_abs(arb_radref(c), arb_midref(a), arb_radref(b));
                _fmpz_addmul_abs(arb_radref(c), arb_midref(b), arb_radref(a));
            }
            else
            {
                /* aliasing */
                fmpz_t t;
                fmpz_init(t);
                fmpz_mul(t, arb_radref(a), arb_radref(b));
                _fmpz_addmul_abs(t, arb_midref(a), arb_radref(b));
                _fmpz_addmul_abs(t, arb_midref(b), arb_radref(a));
                fmpz_swap(arb_radref(c), t);
                fmpz_clear(t);
            }
        }
    }

    fmpz_mul(arb_midref(c), arb_midref(a), arb_midref(b));

    /* common case (both integers) */
    if (fmpz_is_zero(arb_expref(a)) && fmpz_is_zero(arb_expref(b)))
        fmpz_zero(arb_expref(c));
    else
        fmpz_add(arb_expref(c), arb_expref(a), arb_expref(b));

    _arb_normalise(c);
}
