/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sqr(arb_t res, const arb_t x, slong prec)
{
    mag_t resr, xm;
    int inexact;

    if (arb_is_exact(x))
    {
        inexact = arf_sqr_rnd_down(arb_midref(res), arb_midref(x), prec);

        if (inexact)
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(res))
    {
        mag_fast_init_set_arf(xm, arb_midref(x));

        mag_init(resr);
        mag_fast_mul(resr, xm, arb_radref(x));
        if (resr->exp < COEFF_MAX && resr->exp >= COEFF_MIN)
            (resr->exp)++;
        else
            fmpz_add_ui(&(resr->exp), &(resr->exp), 1);
        mag_fast_addmul(resr, arb_radref(x), arb_radref(x));

        inexact = arf_sqr_rnd_down(arb_midref(res), arb_midref(x), prec);

        if (inexact)
            arf_mag_fast_add_ulp(resr, resr, arb_midref(res), prec);

        *arb_radref(res) = *resr;
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));

        mag_init(resr);
        mag_mul(resr, xm, arb_radref(x));
        fmpz_add_ui(&(resr->exp), &(resr->exp), 1);
        mag_addmul(resr, arb_radref(x), arb_radref(x));

        inexact = arf_sqr_rnd_down(arb_midref(res), arb_midref(x), prec);

        if (inexact)
            arf_mag_add_ulp(arb_radref(res), resr, arb_midref(res), prec);
        else
            mag_swap(arb_radref(res), resr);

        mag_clear(xm);
        mag_clear(resr);
    }
}
