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
_fmpz_addmul_abs_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c1 = *g;

    /* c1 is small */
    if (!COEFF_IS_MPZ(c1))
    {
        mp_limb_t prod[2];
        ulong uc1 = FLINT_ABS(c1);

        /* compute product */
        umul_ppmm(prod[1], prod[0], uc1, x);

        /* product fits in one limb */
        if (prod[1] == 0)
        {
            fmpz_add_ui(f, f, prod[0]);
        }
        else
        {
            __mpz_struct * mpz_ptr = _fmpz_promote_val(f);
            mpz_t temp;  /* set up a temporary, cheap mpz_t to contain prod */
            temp->_mp_d = prod;
            temp->_mp_size = 2;
            mpz_add(mpz_ptr, mpz_ptr, temp);
        }
    }
    else  /* c1 is large */
    {
        __mpz_struct * cptr = COEFF_TO_PTR(c1);
        __mpz_struct * mpz_ptr = _fmpz_promote_val(f);
        if (mpz_sgn(cptr) == -1)
            mpz_submul_ui(mpz_ptr, cptr, x);
        else
            mpz_addmul_ui(mpz_ptr, cptr, x);
    }
}

void
_fmpz_addmul_abs(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1, c2;
    __mpz_struct * mpz_ptr, *p1, *p2;

    c1 = *g;
    if (!COEFF_IS_MPZ(c1))  /* g is small */
    {
        _fmpz_addmul_abs_ui(f, h, FLINT_ABS(c1));
        return;
    }

    c2 = *h;
    if (!COEFF_IS_MPZ(c2))  /* h is small */
    {
        _fmpz_addmul_abs_ui(f, g, FLINT_ABS(c2));
        return;
    }

    /* both g and h are large */
    mpz_ptr = _fmpz_promote_val(f);

    p1 = COEFF_TO_PTR(c1);
    p2 = COEFF_TO_PTR(c2);

    if (mpz_sgn(p1) != mpz_sgn(p2))
        mpz_submul(mpz_ptr, p1, p2);
    else
        mpz_addmul(mpz_ptr, p1, p2);
}