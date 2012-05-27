#include "mprb.h"

void
mprb_clear(mprb_t x)
{
    free(x->d);
}