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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "dlog.h"
#include "padic.h"

ulong
dlog_modpe(const dlog_modpe_t t, ulong b)
{
    ulong x;
    x = dlog_precomp(t->modp, b % t->p);
    /*b = b * n_powmod(t->a, -x, t->pe);*/
    b = n_powmod(b, t->p - 1, t->pe);
    x = x + (t->p-1) * dlog_1modpe(t->modpe, b);
    return x;
}
