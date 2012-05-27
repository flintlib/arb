#include "mprb.h"

void
mprb_randtest(mprb_t x, flint_rand_t state, long emin, long emax)
{
    long n;

    n = n_randint(state, x->alloc) + 1;

    _mpr_randtest(x->d, state, n);

    x->rad = n_randtest(state);
    x->size = n;
    x->sign = n_randint(state, 2) ? MPRB_SIGN_PLUS : MPRB_SIGN_MINUS;
    x->exp = emin + n_randint(state, emax - emin + 1);
}
