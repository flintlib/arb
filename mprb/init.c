#include "mprb.h"

void
mprb_init(mprb_t x, long bits)
{
    long limbs = _MPR_BITS_TO_LIMBS(bits);

    x->d = calloc(limbs, sizeof(mp_limb_t));
    x->exp = 0;
    x->sign = MPRB_SIGN_PLUS;
    x->alloc = limbs;
    x->size = limbs;
    x->rad = 0UL;
}
