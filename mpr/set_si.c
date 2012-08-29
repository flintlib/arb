#include "mpr.h"

void mpr_set_si(mpr_t y, long x)
{
    if (x < 0)
    {
        mpr_set_ui(y, -(ulong) x);
        y->sign = MPR_SIGN_MINUS;
    }
    else
    {
        mpr_set_ui(y, x);
    }
}
