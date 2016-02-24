/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* todo: modular arithmetic */

long
n_dirichlet_chi_conrey(const acb_dirichlet_group_t G, const acb_conrey_t a, const acb_conrey_t b)
{
  ulong x, k;
  x = 0;
  for (k = 0; k < G->num; k++)
    x = (x + G->PHI[k] * a->log[k] * b->log[k]) % G->expo;
  return x;
}

long
n_dirichlet_chi(const acb_dirichlet_group_t G, ulong m, ulong n)
{
  ulong x;
  acb_conrey_t a, b;
  acb_conrey_init(a, G);
  acb_conrey_init(b, G);

  acb_conrey_log(a, G, m);
  acb_conrey_log(b, G, n);
  x = n_dirichlet_chi_conrey(G, a, b);

  acb_conrey_clear(a);
  acb_conrey_clear(b);

  return x;
}

void
acb_dirichlet_chi_conrey(acb_t res, const acb_dirichlet_group_t G, const acb_conrey_t a, const acb_conrey_t b, slong prec)
{
  fmpq_t t;
  fmpq_init(t);
  fmpq_set_si(t, n_dirichlet_chi_conrey(G, a, b), G->expo);
  arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), t, prec);
  fmpq_clear(t);
}

void
acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec)
{
  fmpq_t t;
  fmpq_init(t);
  fmpq_set_si(t, n_dirichlet_chi(G, m, n), G->expo);
  arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), t, prec);
  fmpq_clear(t);
}
