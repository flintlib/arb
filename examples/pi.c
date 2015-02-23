/* This file is public domain. Author: Fredrik Johansson. */

#include "arb.h"
#include "profiler.h"

int main(int argc, char *argv[])
{
    arb_t x;
    long prec, digits, condense;

    if (argc < 2)
    {
        printf("usage: build/examples/pi digits [condense = 20]\n");
        return 1;
    }

    digits = atol(argv[1]);

    if (argc > 2)
        condense = atol(argv[2]);
    else
        condense = 20;

    arb_init(x);

    prec = digits * 3.3219280948873623 + 5;

    printf("computing pi with a precision of %ld bits... ", prec);

    TIMEIT_ONCE_START
    arb_const_pi(x, prec);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    arb_printn(x, digits, ARB_STR_CONDENSE * condense);

    printf("\n");

    arb_clear(x);
    flint_cleanup();
    return 0;
}

