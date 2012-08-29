#include "mprb.h"

void
mprb_init(mprb_t x, long bits)
{
    mpr_init(&x->mid);
    ufloat_zero(&x->rad);
    x->bits = bits;
}
