/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <ctype.h>
#include "acb.h"
#include "acb_poly.h"
#include "flint/arith.h"
#include "flint/profiler.h"

int check_accuracy(acb_ptr vec, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (mag_cmp_2exp_si(arb_radref(acb_realref(vec + i)), -prec) >= 0
         || mag_cmp_2exp_si(arb_radref(acb_imagref(vec + i)), -prec) >= 0)
            return 0;
    }

    return 1;
}

slong fmpz_poly_deflation(const fmpz_poly_t input)
{
    slong i, coeff, deflation;

    if (input->length <= 1)
        return input->length;

    coeff = 1;
    while (fmpz_is_zero(input->coeffs + coeff))
        coeff++;

    deflation = n_gcd(input->length - 1, coeff);

    while ((deflation > 1) && (coeff + deflation < input->length))
    {
        for (i = 0; i < deflation - 1; i++)
        {
            coeff++;
            if (!fmpz_is_zero(input->coeffs + coeff))
                deflation = n_gcd(coeff, deflation);
        }

        if (i == deflation - 1)
            coeff++;
    }

    return deflation;
}

void
fmpz_poly_deflate(fmpz_poly_t result, const fmpz_poly_t input, ulong deflation)
{
    slong res_length, i;

    if (deflation == 0)
    {
        flint_printf("Exception (fmpz_poly_deflate). Division by zero.\n");
        abort();
    }

    if (input->length <= 1 || deflation == 1)
    {
        fmpz_poly_set(result, input);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    fmpz_poly_fit_length(result, res_length);
    for (i = 0; i < res_length; i++)
        fmpz_set(result->coeffs + i, input->coeffs + i*deflation);

    result->length = res_length;
}

void
fmpz_poly_complex_roots_squarefree(const fmpz_poly_t poly,
    slong initial_prec,
    slong target_prec,
    slong print_digits)
{
    slong i, j, prec, deg, deg_deflated, isolated, maxiter, deflation;
    acb_poly_t cpoly, cpoly_deflated;
    fmpz_poly_t poly_deflated;
    acb_ptr roots, roots_deflated;
    int removed_zero;

    if (fmpz_poly_degree(poly) < 1)
        return;

    fmpz_poly_init(poly_deflated);
    acb_poly_init(cpoly);
    acb_poly_init(cpoly_deflated);

    /* try to write poly as poly_deflated(x^deflation), possibly multiplied by x */
    removed_zero = fmpz_is_zero(poly->coeffs);
    if (removed_zero)
        fmpz_poly_shift_right(poly_deflated, poly, 1);
    else
        fmpz_poly_set(poly_deflated, poly);
    deflation = fmpz_poly_deflation(poly_deflated);
    fmpz_poly_deflate(poly_deflated, poly_deflated, deflation);

    deg = fmpz_poly_degree(poly);
    deg_deflated = fmpz_poly_degree(poly_deflated);

    flint_printf("searching for %wd roots, %wd deflated\n", deg, deg_deflated);

    roots = _acb_vec_init(deg);
    roots_deflated = _acb_vec_init(deg_deflated);

    for (prec = initial_prec; ; prec *= 2)
    {
        acb_poly_set_fmpz_poly(cpoly_deflated, poly_deflated, prec);
        maxiter = FLINT_MIN(FLINT_MAX(deg_deflated, 32), prec);

        TIMEIT_ONCE_START
        flint_printf("prec=%wd: ", prec);
        isolated = acb_poly_find_roots(roots_deflated, cpoly_deflated,
            prec == initial_prec ? NULL : roots_deflated, maxiter, prec);
        flint_printf("%wd isolated roots | ", isolated);
        TIMEIT_ONCE_STOP

        if (isolated == deg_deflated)
        {
            if (!check_accuracy(roots_deflated, deg_deflated, target_prec))
                continue;

            if (deflation == 1)
            {
                _acb_vec_set(roots, roots_deflated, deg_deflated);
            }
            else  /* compute all nth roots */
            {
                acb_t w, w2;

                acb_init(w);
                acb_init(w2);

                acb_unit_root(w, deflation, prec);
                acb_unit_root(w2, 2 * deflation, prec);

                for (i = 0; i < deg_deflated; i++)
                {
                    if (arf_sgn(arb_midref(acb_realref(roots_deflated + i))) > 0)
                    {
                        acb_root_ui(roots + i * deflation,
                                    roots_deflated + i, deflation, prec);
                    }
                    else
                    {
                        acb_neg(roots + i * deflation, roots_deflated + i);
                        acb_root_ui(roots + i * deflation,
                            roots + i * deflation, deflation, prec);
                        acb_mul(roots + i * deflation,
                            roots + i * deflation, w2, prec);
                    }

                    for (j = 1; j < deflation; j++)
                    {
                        acb_mul(roots + i * deflation + j,
                                roots + i * deflation + j - 1, w, prec);
                    }
                }

                acb_clear(w);
                acb_clear(w2);
            }

            /* by assumption that poly is squarefree, must be just one */
            if (removed_zero)
                acb_zero(roots + deg_deflated * deflation);

            if (!check_accuracy(roots, deg, target_prec))
                continue;

            acb_poly_set_fmpz_poly(cpoly, poly, prec);

            if (!acb_poly_validate_real_roots(roots, cpoly, prec))
                continue;

            for (i = 0; i < deg; i++)
            {
                if (arb_contains_zero(acb_imagref(roots + i)))
                    arb_zero(acb_imagref(roots + i));
            }

            flint_printf("done!\n");
            break;
        }
    }

    if (print_digits != 0)
    {
        _acb_vec_sort_pretty(roots, deg);

        for (i = 0; i < deg; i++)
        {
            acb_printn(roots + i, print_digits, 0);
            flint_printf("\n");
        }
    }

    fmpz_poly_clear(poly_deflated);
    acb_poly_clear(cpoly);
    acb_poly_clear(cpoly_deflated);
    _acb_vec_clear(roots, deg);
    _acb_vec_clear(roots_deflated, deg_deflated);
}

int main(int argc, char *argv[])
{
    fmpz_poly_t f, g;
    fmpz_poly_factor_t fac;
    fmpz_t t;
    slong compd, printd, i, j;

    if (argc < 2)
    {
        flint_printf("poly_roots [-refine d] [-print d] <poly>\n\n");

        flint_printf("Isolates all the complex roots of a polynomial with integer coefficients.\n\n");

        flint_printf("If -refine d is passed, the roots are refined to an absolute tolerance\n");
        flint_printf("better than 10^(-d). By default, the roots are only computed to sufficient\n");
        flint_printf("accuracy to isolate them. The refinement is not currently done efficiently.\n\n");

        flint_printf("If -print d is passed, the computed roots are printed to d decimals.\n");
        flint_printf("By default, the roots are not printed.\n\n");

        flint_printf("The polynomial can be specified by passing the following as <poly>:\n\n");

        flint_printf("a <n>          Easy polynomial 1 + 2x + ... + (n+1)x^n\n");
        flint_printf("t <n>          Chebyshev polynomial T_n\n");
        flint_printf("u <n>          Chebyshev polynomial U_n\n");
        flint_printf("p <n>          Legendre polynomial P_n\n");
        flint_printf("c <n>          Cyclotomic polynomial Phi_n\n");
        flint_printf("s <n>          Swinnerton-Dyer polynomial S_n\n");
        flint_printf("b <n>          Bernoulli polynomial B_n\n");
        flint_printf("w <n>          Wilkinson polynomial W_n\n");
        flint_printf("e <n>          Taylor series of exp(x) truncated to degree n\n");
        flint_printf("m <n> <m>      The Mignotte-like polynomial x^n + (100x+1)^m, n > m\n");
        flint_printf("coeffs <c0 c1 ... cn>        c0 + c1 x + ... + cn x^n\n\n");

        flint_printf("Concatenate to multiply polynomials, e.g.: p 5 t 6 coeffs 1 2 3\n");
        flint_printf("for P_5(x)*T_6(x)*(1+2x+3x^2)\n\n");

        return 1;
    }

    compd = 0;
    printd = 0;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_init(t);
    fmpz_poly_one(f);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-refine"))
        {
            compd = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-print"))
        {
            printd = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "a"))
        {
            slong n = atol(argv[i+1]);
            fmpz_poly_zero(g);
            for (j = 0; j <= n; j++)
                fmpz_poly_set_coeff_ui(g, j, j+1);
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "t"))
        {
            arith_chebyshev_t_polynomial(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "u"))
        {
            arith_chebyshev_u_polynomial(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "p"))
        {
            fmpq_poly_t h;
            fmpq_poly_init(h);
            arith_legendre_polynomial(h, atol(argv[i+1]));
            fmpq_poly_get_numerator(g, h);
            fmpz_poly_mul(f, f, g);
            fmpq_poly_clear(h);
            i++;
        }
        else if (!strcmp(argv[i], "c"))
        {
            arith_cyclotomic_polynomial(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "s"))
        {
            arith_swinnerton_dyer_polynomial(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "b"))
        {
            fmpq_poly_t h;
            fmpq_poly_init(h);
            arith_bernoulli_polynomial(h, atol(argv[i+1]));
            fmpq_poly_get_numerator(g, h);
            fmpz_poly_mul(f, f, g);
            fmpq_poly_clear(h);
            i++;
        }
        else if (!strcmp(argv[i], "w"))
        {
            slong n = atol(argv[i+1]);
            fmpz_poly_zero(g);
            fmpz_poly_fit_length(g, n+2);
            arith_stirling_number_1_vec(g->coeffs, n+1, n+2);
            _fmpz_poly_set_length(g, n+2);
            fmpz_poly_shift_right(g, g, 1);
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "e"))
        {
            fmpq_poly_t h;
            fmpq_poly_init(h);
            fmpq_poly_set_coeff_si(h, 0, 0);
            fmpq_poly_set_coeff_si(h, 1, 1);
            fmpq_poly_exp_series(h, h, atol(argv[i+1]) + 1);
            fmpq_poly_get_numerator(g, h);
            fmpz_poly_mul(f, f, g);
            fmpq_poly_clear(h);
            i++;
        }
        else if (!strcmp(argv[i], "m"))
        {
            fmpz_poly_zero(g);
            fmpz_poly_set_coeff_ui(g, 0, 1);
            fmpz_poly_set_coeff_ui(g, 1, 100);
            fmpz_poly_pow(g, g,  atol(argv[i+2]));
            fmpz_poly_set_coeff_ui(g, atol(argv[i+1]), 1);
            fmpz_poly_mul(f, f, g);
            i += 2;
        }
        else if (!strcmp(argv[i], "coeffs"))
        {
            fmpz_poly_zero(g);
            i++;
            j = 0;
            while (i < argc)
            {
                if (fmpz_set_str(t, argv[i], 10) != 0)
                {
                    i--;
                    break;
                }

                fmpz_poly_set_coeff_fmpz(g, j, t);
                i++;
                j++;
            }
            fmpz_poly_mul(f, f, g);
        }
    }

    fmpz_poly_factor_init(fac);

    flint_printf("computing squarefree factorization...\n");
    TIMEIT_ONCE_START
    fmpz_poly_factor_squarefree(fac, f);
    TIMEIT_ONCE_STOP

    TIMEIT_ONCE_START
    for (i = 0; i < fac->num; i++)
    {
        flint_printf("roots with multiplicity %wd\n", fac->exp[i]);
        fmpz_poly_complex_roots_squarefree(fac->p + i,
            32, compd * 3.32193 + 2, printd);
    }
    TIMEIT_ONCE_STOP

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_clear(t);

    flint_cleanup();
    return EXIT_SUCCESS;
}

