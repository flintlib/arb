/*
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2016 Ricky E. Farr

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_falling(acb_t y, const acb_t x, const acb_t n, slong prec)
{
    if (acb_is_int(n) && arf_sgn(arb_midref(acb_realref(n))) >= 0 &&
        arf_cmpabs_ui(arb_midref(acb_realref(n)), FLINT_MAX(prec, 100)) < 0)
    {
        acb_falling_ui_rec(y, x,
            arf_get_si(arb_midref(acb_realref(n)), ARF_RND_DOWN), prec);
    }
    else
    {
        /* y = gamma(x + 1)/gamma(x - n + 1)  */
        acb_t tmp;
        acb_init(tmp);
        acb_add_ui(tmp, x, 1, prec);
        acb_gamma(y, tmp, prec);
        acb_sub(tmp, tmp, n, prec);
        acb_rgamma(tmp, tmp, prec);
        acb_mul(y, y, tmp, prec);
        acb_clear(tmp);
    }
}
