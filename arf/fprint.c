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

    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "arf.h"

void
arf_fprint(FILE * file, const arf_t x)
{
    if (arf_is_normal(x))
    {
        fmpz_t man, exp;

        fmpz_init(man);
        fmpz_init(exp);

        arf_get_fmpz_2exp(man, exp, x);

        flint_fprintf(file, "(");
        fmpz_fprint(file, man);
        flint_fprintf(file, " * 2^");
        fmpz_fprint(file, exp);
        flint_fprintf(file, ")");

        fmpz_clear(man);
        fmpz_clear(exp);
    }
    else
    {
        if (arf_is_zero(x)) flint_fprintf(file, "(0)");
        else if (arf_is_pos_inf(x)) flint_fprintf(file, "(+inf)");
        else if (arf_is_neg_inf(x)) flint_fprintf(file, "(-inf)");
        else flint_fprintf(file, "(nan)");
    }
}

void
arf_fprintd(FILE * file, const arf_t x, slong d)
{
    if (arf_is_finite(x) && (ARF_EXP(x) <= MPFR_EMIN_MIN + 1 ||
                             ARF_EXP(x) >= MPFR_EMAX_MAX - 1))
    {
        arf_fprint(file, x);
    }
    else
    {
        mpfr_t t;
        mpfr_init2(t, d * 3.33 + 10);
        mpfr_set_emin(MPFR_EMIN_MIN);
        mpfr_set_emax(MPFR_EMAX_MAX);
        arf_get_mpfr(t, x, MPFR_RNDN);
        mpfr_fprintf(file, "%.*Rg", FLINT_MAX(d, 1), t);
        mpfr_clear(t);
    }
}

