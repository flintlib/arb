/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "arb_calc.h"
#include "profiler.h"

long eval_count = 0;

int
z_function(arb_ptr out, const arb_t inp, void * params, long order, long prec)
{
    arb_struct x[2];

    arb_init(x);
    arb_init(x + 1);

    arb_set(x, inp);
    arb_one(x + 1);

    _arb_poly_riemann_siegel_z_series(out, x, FLINT_MIN(2, order), order, prec);

    arb_clear(x);
    arb_clear(x + 1);

    eval_count++;
    return 0;
}

int
sin_x(arb_ptr out, const arb_t inp, void * params, long order, long prec)
{
    int xlen = FLINT_MIN(2, order);

    arb_set(out, inp);
    if (xlen > 1)
        arb_one(out + 1);

    _arb_poly_sin_series(out, out, xlen, order, prec);

    eval_count++;
    return 0;
}

int
sin_x2(arb_ptr out, const arb_t inp, void * params, long order, long prec)
{
    arb_ptr x;

    int xlen = FLINT_MIN(2, order);
    int ylen = FLINT_MIN(3, order);

    x = _arb_vec_init(xlen);

    arb_set(x, inp);
    if (xlen > 1)
        arb_one(x + 1);

    _arb_poly_mullow(out, x, xlen, x, xlen, ylen, prec);
    _arb_poly_sin_series(out, out, ylen, order, prec);

    _arb_vec_clear(x, xlen);

    eval_count++;
    return 0;
}

int
sin_1x(arb_ptr out, const arb_t inp, void * params, long order, long prec)
{
    arb_ptr x;
    int xlen = FLINT_MIN(2, order);

    x = _arb_vec_init(xlen);

    arb_set(x, inp);
    if (xlen > 1)
        arb_one(x + 1);

    _arb_poly_inv_series(out, x, xlen, order, prec);
    _arb_poly_sin_series(out, out, order, order, prec);

    _arb_vec_clear(x, xlen);

    eval_count++;
    return 0;
}

int main(int argc, char *argv[])
{
    arf_interval_ptr blocks;
    arb_calc_func_t function;
    int * info;
    long digits, low_prec, high_prec, i, num, found_roots, found_unknown;
    long maxdepth, maxeval, maxfound;
    int refine;
    double a, b;
    arf_t C;
    arf_interval_t t, interval;
    arb_t v, w, z;

    if (argc < 4)
    {
        printf("real_roots function a b [-refine d] [-verbose] "
            "[-maxdepth n] [-maxeval n] [-maxfound n] [-prec n]\n");
        printf("available functions:\n");
        printf("  0  Z(x), Riemann-Siegel Z-function\n");
        printf("  1  sin(x)\n");
        printf("  2  sin(x^2)\n");
        printf("  3  sin(1/x)\n");
        return 1;
    }

    switch (atoi(argv[1]))
    {
        case 0:
            function = z_function;
            break;
        case 1:
            function = sin_x;
            break;
        case 2:
            function = sin_x2;
            break;
        case 3:
            function = sin_1x;
            break;
        default:
            printf("require a function 0-3\n");
            return 1;
    }

    a = atof(argv[2]);
    b = atof(argv[3]);

    if (a >= b)
    {
        printf("require a < b!\n");
        return 1;
    }

    refine = 0;
    digits = 0;
    maxdepth = 30;
    maxeval = 100000;
    maxfound = 100000;
    low_prec = 30;

    for (i = 4; i < argc; i++)
    {
        if (!strcmp(argv[i], "-refine"))
        {
            refine = 1;
            digits = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-verbose"))
        {
            arb_calc_verbose = 1;
        }
        else if (!strcmp(argv[i], "-maxdepth"))
        {
            maxdepth = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-maxeval"))
        {
            maxeval = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-maxfound"))
        {
            maxfound = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-prec"))
        {
            low_prec = atol(argv[i+1]);
        }
    }

    high_prec = digits * 3.32192809488736 + 10;
    found_roots = 0;
    found_unknown = 0;

    arf_init(C);
    arf_interval_init(t);
    arf_interval_init(interval);
    arb_init(v);
    arb_init(w);
    arb_init(z);

    arf_set_d(&interval->a, a);
    arf_set_d(&interval->b, b);

    printf("interval: "); arf_interval_printd(interval, 15); printf("\n");
    printf("maxdepth = %ld, maxeval = %ld, maxfound = %ld, low_prec = %ld\n",
        maxdepth, maxeval, maxfound, low_prec);

    TIMEIT_ONCE_START

    num = arb_calc_isolate_roots(&blocks, &info, function,
        NULL, interval, maxdepth, maxeval, maxfound, low_prec);

    for (i = 0; i < num; i++)
    {
        if (info[i] != 1)
        {
            if (arb_calc_verbose)
            {
                printf("unable to count roots in ");
                arf_interval_printd(blocks + i, 15);
                printf("\n");
            }
            found_unknown++;
            continue;
        }

        found_roots++;

        if (!refine)
            continue;

        if (arb_calc_refine_root_bisect(t,
            function, NULL, blocks + i, 5, low_prec)
            != ARB_CALC_SUCCESS)
        {
            printf("warning: some bisection steps failed!\n");
        }

        if (arb_calc_verbose)
        {
            printf("after bisection 1: ");
            arf_interval_printd(t, 15);
            printf("\n");
        }

        if (arb_calc_refine_root_bisect(blocks + i,
            function, NULL, t, 5, low_prec)
            != ARB_CALC_SUCCESS)
        {
            printf("warning: some bisection steps failed!\n");
        }

        if (arb_calc_verbose)
        {
            printf("after bisection 2: ");
            arf_interval_printd(blocks + i, 15);
            printf("\n");
        }

        arf_interval_get_arb(v, t, high_prec);
        arb_calc_newton_conv_factor(C, function, NULL, v, low_prec);

        arf_interval_get_arb(w, blocks + i, high_prec);
        if (arb_calc_refine_root_newton(z, function, NULL,
            w, v, C, 10, high_prec) != ARB_CALC_SUCCESS)
        {
            printf("warning: some newton steps failed!\n");
        }

        printf("refined root:\n");
        arb_printd(z, digits + 2);
        printf("\n\n");
    }

    printf("---------------------------------------------------------------\n");
    printf("Found roots: %ld\n", found_roots);
    printf("Subintervals possibly containing undetected roots: %ld\n", found_unknown);
    printf("Function evaluations: %ld\n", eval_count);

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    for (i = 0; i < num; i++)
        arf_interval_clear(blocks + i);
    flint_free(blocks);
    flint_free(info);

    arf_interval_clear(t);
    arf_interval_clear(interval);
    arf_clear(C);
    arb_clear(v);
    arb_clear(w);
    arb_clear(z);
    flint_cleanup();
    return 0;
}

