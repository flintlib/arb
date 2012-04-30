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

#include <alloca.h>
#include "mpr.h"

/*
Evaluates y = c0 + c1*x + ... + c[len-1]*x^(len-1) where

c0 are limb-size integers
x is a fixed-point number in [0,1) with <prec> limbs
y is a fixed-point number in [0,LIMB_MAX) with <prec>+1 limbs

Assumes no aliasing of x, y, coeffs.
Assumes precisions small enough to allocate everything on the stack.

Temporary allocation:

powers: x^2, x^3, ..., x^m, each of size prec
T: temporary sum of size prec + 1
U: temporary product of size 2*prec + 1

|  T  |  U  |  x^m  |  x^(m-1)  |  ...  |  x^2  |

total size: (prec + 1) + (2*prec + 1) + (m-1)*prec = 2 + (2+m)*prec

With the reverse order of powers and padding to the left of x^m,
we can perform multiplications in-place, leaving the low-order limbs
to be overwritten.

TODO: error analysis.
*/

/* theoretically optimal parameters (minimize number of multiplications,
   and on a tie, minimize the number of multiplications in the second step) */
const unsigned char opt_split[] =
{
    0, 0, 2, 3, 2, 3, 3, 4, 4, 3, 5, 4, 4, 5, 5, 5, 4,
    6, 6, 5, 5, 7, 6, 6, 6, 5, 7, 7, 7, 6, 6, 8, 8, 7,
    7, 7, 6, 8, 8, 8, 8, 7, 7, 9, 9, 9, 8, 8, 8, 7, 10,
    9, 9, 9, 9, 8, 8, 10, 10, 10, 10, 9, 9, 9, 8, 11, 11,
    10, 10, 10, 10, 9, 9, 11, 11, 11, 11, 11, 10, 10, 10,
    9, 12, 12, 12, 11, 11, 11, 11, 10, 10, 13, 12, 12, 12,
    12, 12, 11, 11, 11, 10, 13, 13, 13, 13, 12, 12, 12, 12,
    11, 11, 14, 14, 13, 13, 13, 13, 13, 12, 12, 12, 11, 14,
    14, 14, 14, 14, 13, 13, 13, 13, 12, 12, 15, 15, 15, 14,
    14, 14, 14, 14, 13, 13, 13, 12, 15, 15, 15, 15, 15, 15,
    14, 14, 14, 14, 13, 13, 16, 16, 16, 16, 15, 15, 15, 15,
    15, 14, 14, 14, 13, 17, 16, 16, 16, 16, 16, 16, 15, 15,
    15, 15, 14, 14, 17, 17, 17, 17, 17, 16, 16, 16, 16, 16,
    15, 15, 15, 14, 18, 18, 17, 17, 17, 17, 17, 17, 16, 16,
    16, 16, 15, 15, 18, 18, 18, 18, 18, 18, 17, 17, 17, 17,
    17, 16, 16, 16, 15, 19, 19, 19, 18, 18, 18, 18, 18, 18,
    17, 17, 17, 17, 16, 16, 19, 19, 19, 19, 19, 19, 19, 18,
    18, 18, 18, 18, 17, 17, 17
};

void
mpr_polyval_1(mp_ptr y, mp_srcptr x, long prec, mp_srcptr coeffs, long len)
{
    mp_ptr tmp, T, U, P;
    long i, alloc, m, n1, n2;

    if (len <= 3)
    {
        if (len == 0)
        {
            mpn_zero(y, prec + 1);
        }
        else if (len == 1)
        {
            mpn_zero(y, prec);
            y[prec] = coeffs[0];
        }
        else if (len == 2)
        {
            y[prec] = mpn_mul_1(y, x, prec, coeffs[1]) + coeffs[0];
        }
        else if (len == 3)
        {
            tmp = alloca(2 * prec);
            mpn_sqr(tmp, x, prec);
            y[prec] = mpn_mul_1(y, x, prec, coeffs[1]);
            y[prec] += mpn_addmul_1(y, tmp, prec, coeffs[2]);
            y[prec] += coeffs[0];
        }
        return;
    }

    /* rectangular splitting parameter */
    m = opt_split[len];

    alloc = 2 + (2+m)*prec;
    tmp = alloca(alloc * sizeof(mp_limb_t));
    T = tmp;
    U = T + prec + 1;
    P = U + (2*prec + 1);

#define POW_WRITE(j) (P + (m - (j)) * prec - prec)
#define POW_READ(j)  (P + (m - (j)) * prec)

    /* compute powers */
    mpn_sqr(POW_WRITE(2), x, prec);
    for (i = 3; i <= m; i++)
        mpn_mul_n(POW_WRITE(i), POW_READ(i-1), x, prec);

    /* evaluate the first row */
    n2 = len - 1;
    n1 = m * (n2 / m);
    if (n2 == n1)
    {
        mpn_zero(y, prec);
        y[prec] = coeffs[n1];
    }
    else
    {
        y[prec] = mpn_mul_1(y, x, prec, coeffs[n1 + 1]) + coeffs[n1];
        for (i = n1 + 2; i <= n2; i++)
            y[prec] += mpn_addmul_1(y, POW_READ(i - n1), prec, coeffs[i]);
    }

    /* evaluate remaining rows */
    while (n2 > m)
    {
        n2 = n1;
        n1 = n2 - m;
        n2 -= 1;

        T[prec] = mpn_mul_1(T, x, prec, coeffs[n1 + 1]) + coeffs[n1];
        for (i = n1 + 2; i <= n2; i++)
            T[prec] += mpn_addmul_1(T, POW_READ(i - n1), prec, coeffs[i]);

        /* U = y * x^m */
        mpn_mul(U, y, prec + 1, POW_READ(m), prec);
        /* y = T + U */
        mpn_add_n(y, T, U + prec, prec + 1);
    }
}

