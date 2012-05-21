#include "mpr.h"

void
_mpr_printr(mp_ptr x, mp_size_t n)
{
    long i;

    printf("[");
    for (i = 0; i < n; i++)
    {
        printf("0x%lx", x[i]);
        if (i < n - 1)
            printf(", ");
    }

    printf("]\n");
}
