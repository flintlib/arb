/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "ufloat.h"

/* must be smaller than UFLOAT_BITS */
#define LOG_FRAC_BITS 8
#define LOG_INDEX_BITS 8

#define LOG_INDEX_OFFSET 128
#define LOG_VALUE_OFFSET 1243UL

#define LOG2_LBOUND 177UL
#define LOG2_UBOUND 178UL

/* [int(ceil(log(n) * 2**8))-LOG_MIN for n in range(128,257)] */
unsigned char ufloat_log_tab[] = {
    0, 2, 4, 6, 7, 9, 11, 13, 15, 17, 19, 21, 23, 24, 26, 28, 30, 32,
    33, 35, 37, 39, 40, 42, 44, 45, 47, 49, 50, 52, 54, 55, 57, 58, 60,
    62, 63, 65, 66, 68, 69, 71, 72, 74, 75, 77, 78, 80, 81, 83, 84, 85,
    87, 88, 90, 91, 93, 94, 95, 97, 98, 99, 101, 102, 103, 105, 106, 107,
    109, 110, 111, 113, 114, 115, 116, 118, 119, 120, 121, 123, 124, 125,
    126, 128, 129, 130, 131, 132, 134, 135, 136, 137, 138, 139, 141, 142,
    143, 144, 145, 146, 147, 149, 150, 151, 152, 153, 154, 155, 156, 157,
    158, 159, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
    173, 174, 175, 176, 177
};

void
ufloat_log(ufloat_t z, const ufloat_t x)
{
    ufloat_t t, u;

    /* log(1) = 0 */
    if ((x->exp == (1L - UFLOAT_PREC) && x->man == (1UL << (UFLOAT_PREC - 1)))
        || (x->exp <= -UFLOAT_PREC))
    {
        z->man = 0;
        z->exp = 0;
        return;
    }

    t->man = n_rshift_ceil(x->man, UFLOAT_PREC - LOG_INDEX_BITS);
    t->man = (mp_limb_t) ufloat_log_tab[t->man - LOG_INDEX_OFFSET] + LOG_VALUE_OFFSET;
    t->exp = -LOG_FRAC_BITS;
    ufloat_normalise(t);

    /* log(m * 2^e) = log(m) + e*log(2) */
    if ((x->exp + (UFLOAT_PREC - LOG_INDEX_BITS)) >= 0)
    {
        u->man = LOG2_UBOUND * (x->exp + (UFLOAT_PREC - LOG_INDEX_BITS));
        u->exp = -LOG_FRAC_BITS;
        ufloat_normalise(u);
        ufloat_add(z, t, u);
        return;
    }
    else
    {
        u->man = LOG2_LBOUND * (-(x->exp + (UFLOAT_PREC - LOG_INDEX_BITS)));
        u->exp = -LOG_FRAC_BITS;
        ufloat_normalise(u);
        ufloat_sub(z, t, u);
        return;
    }
}
