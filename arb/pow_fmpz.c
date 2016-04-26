/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, slong prec)
{
    arb_pow_fmpz_binexp(y, b, e, prec);
}

void
arb_pow_ui(arb_t y, const arb_t b, ulong e, slong prec)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    arb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}

void
arb_ui_pow_ui(arb_t y, ulong b, ulong e, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, b);
    arb_pow_ui(y, t, e, prec);
    arb_clear(t);
}

void
arb_si_pow_ui(arb_t y, slong b, ulong e, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_si(t, b);
    arb_pow_ui(y, t, e, prec);
    arb_clear(t);
}

