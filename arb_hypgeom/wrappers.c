/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
arb_hypgeom_erf(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_erf(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_erfc(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_erfc(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_erfi(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_erfi(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, slong prec)
{
    if (!arb_is_finite(z))
    {
        if (res1 != NULL) arb_indeterminate(res1);
        if (res2 != NULL) arb_indeterminate(res2);
    }
    else
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);
        arb_set(acb_realref(t), z);
        acb_hypgeom_fresnel(res1 ? t : NULL, res2 ? u : NULL, t, normalized, prec);
        if (res1 != NULL) arb_swap(res1, acb_realref(t));
        if (res2 != NULL) arb_swap(res2, acb_realref(u));
        acb_clear(t);
        acb_clear(u);
    }
}

