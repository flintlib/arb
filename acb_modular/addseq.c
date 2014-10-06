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

#include "acb_modular.h"

static long
bisect(long needle, const long * haystack, long len)
{
    long a, b, mid;

    a = 0;
    b = len - 1;

    while (a < b)
    {
        mid = (a + b) / 2;

        if (haystack[mid] < needle)
            a = mid + 1;
        else
            b = mid;
    }

    if (a == b && haystack[a] == needle)
    {
        return a;
    }

    return -1;
}

/* write p = 2a with a in P */
static int
write_as_2a(long * i1, long * i2, long p, const long * P, long Plen)
{
    if (p % 2 == 0)
    {
        long i = bisect(p / 2, P, Plen);

        if (i != -1)
        {
            *i1 = *i2 = i;
            return 1;
        }
    }

    return 0;
}

/* write p = a + b with a, b in P */
static int
write_as_a_b(long * i1, long * i2, long p, const long * P, long Plen)
{
    long i, j, pi;

    for (i = 0; i < Plen; i++)
    {
        pi = P[i];

        j = bisect(p - pi, P, Plen);

        if (j != -1)
        {
            *i1 = i;
            *i2 = j;
            return 1;
        }
    }

    return 0;
}

/* write p = 2a + b with a, b in P */
static int
write_as_2a_b(long * i1, long * i2, long p, const long * P, long Plen)
{
    long i, j, pi;

    for (i = 0; i < Plen; i++)
    {
        pi = P[i];

        if (2 * pi > p)
            break;

        j = bisect(p - 2*pi, P, Plen);

        if (j != -1)
        {
            *i1 = i;
            *i2 = j;
            return 1;
        }
    }

    return 0;
}

void
acb_modular_addseq_theta(long * exponents, long * aindex, long * bindex, long num)
{
    long i;
    long c;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
        {
            exponents[i] = 1;
            aindex[i] = 0;
            bindex[i] = 0;
            continue;
        }
        else
        {
            c = ((i + 2) * (i + 2)) / 4;

            exponents[i] = c;

            if (write_as_2a(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_2a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            printf("i = %ld, c = %lu: bad addition sequence!\n", i, c);
            abort();
        }
    }
}

void
acb_modular_addseq_eta(long * exponents, long * aindex, long * bindex, long num)
{
    long i;
    long c;

    for (i = 0; i < num; i++)
    {
        if (i == 0)
        {
            exponents[i] = 1;
            aindex[i] = 0;
            bindex[i] = 0;
            continue;
        }
        else
        {
            c = ((i + 2) / 2) * ((3 * i + 5) / 2) / 2;

            exponents[i] = c;

            if (write_as_2a(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            if (write_as_2a_b(aindex + i, bindex + i, c, exponents, i))
                continue;

            printf("i = %ld, c = %lu: bad addition sequence!\n", i, c);
            abort();
        }
    }
}

