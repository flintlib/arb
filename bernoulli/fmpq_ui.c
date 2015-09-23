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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "bernoulli.h"

void
_bernoulli_fmpq_ui(fmpz_t num, fmpz_t den, ulong n)
{
    if (n < (ulong) bernoulli_cache_num)
    {
        fmpz_set(num, fmpq_numref(bernoulli_cache + n));
        fmpz_set(den, fmpq_denref(bernoulli_cache + n));
    }
    else
    {
        _bernoulli_fmpq_ui_zeta(num, den, n);
    }
}

void
bernoulli_fmpq_ui(fmpq_t b, ulong n)
{
    _bernoulli_fmpq_ui(fmpq_numref(b), fmpq_denref(b), n);
}

