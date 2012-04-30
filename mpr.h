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

#ifndef MPR_H
#define MPR_H

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include <math.h>
#include "flint.h"
#include "ulong_extras.h"

void mpn_mul_basecase(mp_ptr z, mp_srcptr xp, mp_size_t xn, mp_srcptr yp, mp_size_t yn);

void mpr_polyval_1(mp_ptr y, mp_srcptr x, long prec, mp_srcptr coeffs, long len);

void mpr_exp_basecase(mp_ptr y, mp_srcptr x, long prec, mp_bitcnt_t bits);


#endif
