/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define N 10000

int main()
{
  slong iter;
  slong prec;
  flint_rand_t state;
  arf_rnd_t round;
  arb_t rand[N];
  arb_t m; /* mean */
  arb_t s; /* variance */
  arb_t mp;
  arb_t sp;

  flint_printf("urandom....");
  fflush(stdout);

  flint_randinit(state);
  arb_init(m);
  arb_init(s);
  arb_init(mp);
  arb_init(sp);
  prec = 299;
  round = ARF_RND_DOWN;

  for (iter = 0; iter < N; iter++)
  {
    arb_init(rand[iter]);
    arb_urandom(rand[iter], state, prec, round);
    arb_add(m, m, rand[iter], prec);
  }

  arb_div_si(m, m, N, prec);

  for (iter = 0; iter < N; iter++)
  {
    arb_t tmp;
    arb_init(tmp);
    arb_sub(tmp, rand[iter], m, prec);
    arb_sqr(tmp, tmp, prec);
    arb_add(s, s, tmp, prec);
    arb_clear(tmp);
  }

  arb_div_si(s, s, N, prec);

  /* one percent deviation */
  arb_set_str(mp, "0.5 +/- 0.005", prec);
  arb_set_str(sp, "0.083333 +/- 0.00083", prec);

  if (!arb_contains(mp, m))
  {
      flint_printf("FAIL: mean\n\n");
      flint_printf("m = "); arb_printd(m, 15); flint_printf("\n\n");
      flint_abort();
  }
  if (!arb_contains(sp, s))
  {
      flint_printf("FAIL: variance\n\n");
      flint_printf("s = "); arb_printd(s, 15); flint_printf("\n\n");
      flint_abort();
  }

  for (iter = 0; iter < N; iter++) arb_clear(rand[iter]);
  arb_clear(m);
  arb_clear(s);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
