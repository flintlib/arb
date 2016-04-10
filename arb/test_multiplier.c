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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "arb.h"

double _arb_test_multiplier = -1.0;

double arb_test_multiplier()
{
    if (_arb_test_multiplier == -1.0)
    {
        const char * s = getenv("ARB_TEST_MULTIPLIER");

        if (s == NULL)
        {
            _arb_test_multiplier = 1.0;
        }
        else
        {
            _arb_test_multiplier = strtod(s, NULL);

            if (!(_arb_test_multiplier >= 0.0 && _arb_test_multiplier <= 1000.0))
                _arb_test_multiplier = 1.0;
        }
    }

    return _arb_test_multiplier;
}

