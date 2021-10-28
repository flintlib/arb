/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "bernoulli.h"
#include "arb.h"

#define TIMING 0
#define DEBUG 0

static void
_fmpz_crt_combine(fmpz_t r1r2, fmpz_t m1m2, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2)
{
    fmpz_invmod(m1m2, m1, m2);
    fmpz_mul(m1m2, m1m2, m1);
    fmpz_sub(r1r2, r2, r1);
    fmpz_mul(r1r2, r1r2, m1m2);
    fmpz_add(r1r2, r1r2, r1);
    fmpz_mul(m1m2, m1, m2);
    fmpz_mod(r1r2, r1r2, m1m2);
}

static void
tree_crt(fmpz_t r, fmpz_t m, mp_srcptr residues, mp_srcptr primes, slong len)
{
    if (len == 0)
    {
        fmpz_zero(r);
        fmpz_one(m);
    }
    else if (len == 1)
    {
        fmpz_set_ui(r, residues[0]);
        fmpz_set_ui(m, primes[0]);
    }
    else
    {
        fmpz_t r1, m1, r2, m2;

        fmpz_init(r1);
        fmpz_init(m1);
        fmpz_init(r2);
        fmpz_init(m2);

        tree_crt(r1, m1, residues, primes, len / 2);
        tree_crt(r2, m2, residues + len / 2, primes + len / 2, len - len / 2);
        _fmpz_crt_combine(r, m, r1, m1, r2, m2);

        fmpz_clear(r1);
        fmpz_clear(m1);
        fmpz_clear(r2);
        fmpz_clear(m2);
    }
} 

void
_bernoulli_fmpq_ui_multi_mod(fmpz_t num, fmpz_t den, ulong n, double alpha)
{
    slong i, bits, mod_bits, zeta_bits, num_primes;
    ulong p;
    mp_ptr primes, residues;
    mag_t primes_product;
    fmpz_t M;
#if TIMING
    double t1, t2;
#endif

    if (n < 10 || n % 2 != 0)
    {
        _bernoulli_fmpq_ui_zeta(num, den, n);
        return;
    }

    if (alpha < 0)
    {
        if (n < 18000)
            alpha = 0.0;
        else if (n < 60000)
            alpha = 0.005 + 3.6e-6 * n;
        else
            alpha = FLINT_MIN(0.18 + 0.5e-6 * n, 0.28);
    }

#if TIMING
    t1 = clock();
#endif

    arith_bernoulli_number_denom(den, n);

    bits = arith_bernoulli_number_size(n) + fmpz_bits(den) + 2;
    mod_bits = bits * alpha;
    zeta_bits = bits - mod_bits;

    num_primes = 0;
    mag_init(primes_product);
    mag_one(primes_product);

    for (p = 5; mag_cmp_2exp_si(primes_product, mod_bits) < 0; p = n_nextprime(p, 1))
    {
        if (n % (p - 1) != 0)
        {
            mag_mul_ui_lower(primes_product, primes_product, p);
            num_primes++;
        }
    }

#if DEBUG
    printf("\nn = %lu, bits = %lu, num_primes = %ld\n", n, bits, num_primes);
#endif

    primes = flint_malloc(sizeof(mp_limb_t) * num_primes);
    residues = flint_malloc(sizeof(mp_limb_t) * num_primes);

    for (p = 5, i = 0; i < num_primes; p = n_nextprime(p, 1))
    {
        if (n % (p - 1) != 0)
        {
            primes[i] = p;
            i++;
        }
    }

#if TIMING
    t2 = clock();
    printf("init time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    printf("num_primes = %ld\n", num_primes);
#endif

    for (i = 0; i < num_primes; i++)
    {
#if TIMING
        if (i % 10000 == 0)
            printf("%ld / %ld\n", i, num_primes);
#endif

        residues[i] = bernoulli_mod_p_harvey(n, primes[i]);
    }

#if TIMING
    t2 = clock();
    printf("mod time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    printf("start CRT\n");
    t1 = clock();
#endif

    fmpz_init(M);
    tree_crt(num, M, residues, primes, num_primes);
    fmpz_mul(num, num, den);
    fmpz_mod(num, num, M);

    if (n % 4 == 0)
    {
        fmpz_sub(num, M, num);
        fmpz_neg(num, num);
    }

#if TIMING
    printf("end CRT\n");
    t2 = clock();
    printf("CRT time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    t1 = clock();
#endif

    if (zeta_bits > 0)
    {
        slong prec;
        arb_t b;
        fmpz_t t;

        arb_init(b);
        fmpz_init(t);

        for (prec = zeta_bits + 10; ; prec += 32)
        {
            arb_bernoulli_ui_zeta(b, n, prec);
            arb_mul_fmpz(b, b, den, prec);
            arb_sub_fmpz(b, b, num, prec);
            arb_div_fmpz(b, b, M, prec);

            if (arb_get_unique_fmpz(t, b))
            {
                fmpz_addmul(num, t, M);
                break;
            }

            flint_printf("bernoulli: n = %wu, bits = %wd, mod = %wd, zeta = %wd: get_unique_fmpz failed!\n", n, bits, mod_bits, zeta_bits);
        }

        arb_clear(b);
        fmpz_clear(t);
    }

#if TIMING
    printf("end zeta\n");
    t2 = clock();
    printf("zeta time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
#endif

    flint_free(primes);
    flint_free(residues);
    fmpz_clear(M);
    mag_clear(primes_product);
}

