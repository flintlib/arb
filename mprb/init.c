#include "mprb.h"

void
mprb_init(mprb_t x, long bits)
{
    long limbs = _MPR_BITS_TO_LIMBS(bits);

    x->mid.d = calloc(limbs, sizeof(mp_limb_t));
    x->mid.exp = 0;
    x->mid.sign = MPRB_SIGN_PLUS;
    x->mid.alloc = limbs;
    x->mid.size = limbs;

    /* XXX: ufloat_zero */
    x->rad.man = 0UL;
    x->rad.exp = 0L; 
}
