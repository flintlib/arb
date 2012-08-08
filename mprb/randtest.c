#include "mprb.h"

void
mprb_randtest(mprb_t x, flint_rand_t state, long emin, long emax)
{
    long n;

#if 0
    n = x->mid.alloc;
    _mpr_randtest(x->mid.d, state, n);
    x->mid.d[0] |= 1UL;
#else
    n = n_randint(state, x->mid.alloc) + 1;
    _mpr_randtest(x->mid.d, state, n);
#endif

    x->mid.size = n;
    x->mid.sign = n_randint(state, 2) ? MPRB_SIGN_PLUS : MPRB_SIGN_MINUS;
    x->mid.exp = emin + n_randint(state, emax - emin + 1);

    /* XXX: change ufloat_randtest to use emin / emax */
    ufloat_randtest(&x->rad, state, emax);
}
