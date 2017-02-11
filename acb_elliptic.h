/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACB_ELLIPTIC_H
#define ACB_ELLIPTIC_H

#include <stdio.h>
#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec);

void acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec);

void acb_elliptic_rc1(acb_t res, const acb_t x, slong prec);

#ifdef __cplusplus
}
#endif

#endif

