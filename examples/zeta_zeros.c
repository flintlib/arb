/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "acb_dirichlet.h"
#include "flint/profiler.h"

int main(int argc, char *argv[])
{
    arb_ptr res;
    slong i, prec, num, digits, num_threads, found, alloc;
    int platt, suppress, stream;
    fmpz_t start;

    if (argc < 2)
    {
        flint_printf("usage: zeta_zeros [-start n] [-num n] [-threads t] [-digits d] [-suppress] [-stream] [-platt]\n");
        flint_printf("compute nontrivial zeros of the Riemann zeta function\n");
        return 1;
    }

    fmpz_init(start);
    fmpz_one(start);

    digits = 30;
    num_threads = 1;
    suppress = 0;
    num = 1;
    platt = 0;
    stream = 0;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-suppress"))
        {
            suppress = 1;
        }
        if (!strcmp(argv[i], "-stream"))
        {
            stream = 1;
        }
        else if (!strcmp(argv[i], "-num"))
        {
            num = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-start"))
        {
            fmpz_set_str(start, argv[i+1], 10);
            i++;
        }
        else if (!strcmp(argv[i], "-digits"))
        {
            digits = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-platt"))
        {
            platt = 1;
        }
    }

    flint_set_num_threads(num_threads);

    prec = digits * 3.3219280948873623 + 3;

    flint_printf("Computing %wd zeros starting with zero #", num);
    fmpz_print(start);
    flint_printf(", with prec = %wd bits...\n", prec);
    fflush(stdout);

    TIMEIT_ONCE_START
    if (platt)
    {
        alloc = num;
        res = _arb_vec_init(alloc);

        found = acb_dirichlet_platt_hardy_z_zeros(res, start, num, prec);
    }
    else if (stream)
    {
        fmpz_t n;
        fmpz_init(n);

        alloc = 1;
        res = _arb_vec_init(alloc);

        found = 0;
        for (i = 0; i < num; i++)
        {
            fmpz_add_ui(n, start, i);
            acb_dirichlet_hardy_z_zero(res, n, prec);
            if (!suppress)
            {
                arb_printn(res, digits, ARB_STR_NO_RADIUS);
                flint_printf("\n");
            }
            found++;
        }

        fmpz_clear(n);
    }
    else
    {
        alloc = num;
        res = _arb_vec_init(alloc);

        acb_dirichlet_hardy_z_zeros(res, start, num, prec);
        found = num;
    }
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    if (found < num)
        flint_printf("Found %wd zeros\n", found);

    if (!suppress && !stream)
    {
        for (i = 0; i < found; i++)
        {
            arb_printn(res + i, digits, ARB_STR_NO_RADIUS);
            flint_printf("\n");
        }
    }

    flint_printf("\n");

    _arb_vec_clear(res, alloc);
    fmpz_clear(start);
    flint_cleanup_master();
    return 0;
}

