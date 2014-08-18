/* This file is public domain. Author: Fredrik Johansson. */

#include "arb.h"
#include "profiler.h"

int main(int argc, char *argv[])
{
    arb_t x;
    long prec, digits, digits_to_print;

    if (argc < 2)
    {
        printf("usage: build/examples/pi digits [digits_to_print = 10]\n");
        return 1;
    }

    digits = atol(argv[1]);

    if (argc > 2)
        digits_to_print = atol(argv[2]);
    else
        digits_to_print = 10;

    arb_init(x);

    prec = digits * 3.3219280948873623 + 5;

    printf("computing pi with a precision of %ld bits... ", prec);

    TIMEIT_ONCE_START
    arb_const_pi(x, prec);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    arb_printd(x, digits_to_print);

    printf("\n");

    arb_clear(x);
    flint_cleanup();
    return 0;
}

