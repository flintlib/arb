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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "fmpz_extras.h"

void
fmpz_lshift_mpn(fmpz_t z, mp_srcptr d, mp_size_t dn, int sgnbit, mp_bitcnt_t shift)
{
    __mpz_struct * zmpz;
    mp_ptr zp;
    mp_size_t zn, shift_limbs;
    mp_bitcnt_t shift_bits;

    zmpz = _fmpz_promote(z);

    shift_limbs = shift / FLINT_BITS;
    shift_bits = shift % FLINT_BITS;
    zn = dn + shift_limbs + (shift_bits != 0);

    if (zmpz->_mp_alloc < zn)
        mpz_realloc2(zmpz, zn * FLINT_BITS);

    zp = zmpz->_mp_d;
    flint_mpn_zero(zp, shift_limbs);

    if (shift_bits == 0)
    {
        flint_mpn_copyi(zp + shift_limbs, d, dn);
    }
    else
    {
        zp[zn - 1] = mpn_lshift(zp + shift_limbs, d, dn, shift_bits);
        while (zp[zn - 1] == 0)
            zn--;
    }

    zmpz->_mp_size = sgnbit ? -(long) zn : zn;
    _fmpz_demote_val(z);
}

