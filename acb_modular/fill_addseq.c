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

#include "acb_modular.h"

void
acb_modular_fill_addseq(long * tab, long len)
{
    long i, j;

    for (i = 2; i < len; i++)
    {
        if (tab[i] == -1)
        {
            /* prefer doubling (squaring) */
            if ((i % 2) == 0 && tab[i / 2] != 0)
            {
                tab[i] = i / 2;
            }
            else
            {
                /* check if it can be written as a sum */
                /* prefer unbalanced, because one power will have low precision */
#if 1
                for (j = 1; 2 * j <= i; j++)
#else
                for (j = i / 2; j >= 1; j--)
#endif
                {
                    if (tab[j] != 0 && tab[i-j] != 0)
                    {
                        tab[i] = j;
                        break;
                    }
                }

                /* extend and start over */
                if (tab[i] == -1)
                {
                    tab[i] = i / 2;

                    if (tab[i / 2] == 0)
                        tab[i / 2] = -1;

                    if (tab[i - i / 2] == 0)
                        tab[i - i / 2] = -1;

                    i = 1;
                }
            }
        }
    }

    /* prefer squaring (extra entries may have been inserted) */
    for (i = 2; i < len; i += 2)
    {
        if (tab[i] != 0 && tab[i] != i / 2)
        {
            if (tab[i / 2] != 0)
                tab[i] = i / 2;
        }
    }
}

