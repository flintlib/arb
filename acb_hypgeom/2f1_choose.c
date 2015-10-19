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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

#define ALWAYS1 (0.25 * 0.25)
#define ALWAYS2 (0.75 * 0.75)
#define LIMIT (0.75 * 0.75)

int
acb_hypgeom_2f1_choose(const acb_t z)
{
    double x, y;
    double mag[7];
    int i, pick;

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    x = FLINT_MAX(FLINT_MIN(x, 1e10), -1e10);
    y = FLINT_MAX(FLINT_MIN(y, 1e10), -1e10);

    mag[0] = x*x + y*y;  /* |z|^2 */
    mag[4] = (1.0-x)*(1.0-x) + y*y;              /* |1-z|^2 */

    if (mag[0] <= ALWAYS1)   return 0;

    mag[1] = mag[0] / FLINT_MAX(mag[4], 1e-10);  /* |z/(z-1)|^2 */

    if (mag[1] <= ALWAYS1)   return 1;

    if (mag[0] <= ALWAYS2 || mag[1] <= ALWAYS2)
        return mag[0] <= mag[1] ? 0 : 1;

    mag[2] = 1.0 / mag[0];                    /* |1/z|^2 */
    mag[3] = 1.0 / FLINT_MAX(mag[4], 1e-10);  /* 1/|1-z|^2 */
    mag[5] = mag[4] / mag[0];                 /* |1-1/z|^2 = |(1-z)/z|^2 */

    pick = 0;
    for (i = 1; i < 6; i++)
    {
        if (mag[i] < mag[pick])
            pick = i;
    }

    if (mag[pick] <= LIMIT)
        return pick;

    return 6;
}

