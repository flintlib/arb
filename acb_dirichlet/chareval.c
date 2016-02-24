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

    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_dirichlet.h"

long
n_dirichlet_char_eval(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n)
{
  ulong v = 0, k;
  acb_conrey_t x;
  acb_conrey_init(x, G);
  acb_conrey_log(x, G, n);
  for (k = 0; k < G->num; k++)
    v = (v + chi->expo[k] * x->log[k]) % chi->order;
  acb_conrey_clear(x);
  return v;
}

void
fmpq_dirichlet_char_eval(fmpq_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n)
{
    fmpq_set_si(res, n_dirichlet_char_eval(G, chi, n), chi->order);
}

void
acb_dirichlet_char_eval(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n, slong prec)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_dirichlet_char_eval(t, G, chi, n);
    arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), t, prec);
    fmpq_clear(t);
}
