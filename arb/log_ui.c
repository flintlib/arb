/*
    Copyright (C) 2013-2014 Fredrik Johansson
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/thread_support.h"
#include "arb.h"

/*
Hyperbolic arctangent formulas for log(p).
See https://www.jjj.de/fxt/fxtbook.pdf section 32.4.
*/

#define NUM_SMALL_PRIMES 13

/* static const ulong small_primes[NUM_SMALL_PRIMES] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41 }; */

static const ulong small_primes_atanh_x[NUM_SMALL_PRIMES] = { 51744295, 170918749, 265326335,
    287080366, 362074049, 587270881, 617831551, 740512499,
    831409151, 1752438401, 2151548801, 2470954914, 3222617399 };

static const slong small_primes_atanh_c[NUM_SMALL_PRIMES][NUM_SMALL_PRIMES] = {
  { -1595639, -17569128, -8662593, -31112926, -13108464, -11209640, -12907342, +9745611, -1705229, -12058985, +4580610, +4775383, -12972664 },
  { -2529028, -27846409, -13729885, -49312821, -20776424, -17766859, -20457653, +15446428, -2702724, -19113039, +7260095, +7568803, -20561186 },
  { -3704959, -40794252, -20113918, -72241977, -30436911, -26027978, -29969920, +22628608, -3959419, -28000096, +10635847, +11088096, -30121593 },
  { -4479525, -49322778, -24318973, -87345026, -36800111, -31469438, -36235490, +27359389, -4787183, -33853851, +12859398, +13406195, -36418872 },
  { -5520004, -60779197, -29967648, -107633040, -45347835, -38778983, -44652067, +33714275, -5899123, -41717234, +15846307, +16520111, -44878044 },
  { -5904566, -65013499, -32055403, -115131507, -48507081, -41480597, -47762841, +36063046, -6310097, -44623547, +16950271, +17671017, -48004561 },
  { -6522115, -71813158, -35408027, -127172929, -53580360, -45818987, -52758281, +39834823, -6970060, -49290653, +18723073, +19519201, -53025282 },
  { -6778159, -74632382, -36798067, -132165454, -55683805, -47617738, -54829453, +41398649, -7243689, -51225694, +19458099, +20285481, -55106936 },
  { -7217972, -79475039, -39185776, -140741248, -59296949, -50707501, -58387161, +44084875, -7713709, -54549566, +20720673, +21601741, -58682649 },
  { -7751584, -85350490, -42082712, -151146003, -63680669, -54456218, -62703622, +47343993, -8283970, -58582320, +22252516, +23198720, -63020955 },
  { -7905109, -87040909, -42916186, -154139543, -64941904, -55534757, -63945506, +48281670, -8448039, -59742579, +22693241, +23658185, -64269124 },
  { -8312407, -91525553, -45127374, -162081337, -68287932, -58396097, -67240196, +50769306, -8883311, -62820720, +23862474, +24877135, -67580488 },
  { -8548719, -94127517, -46410292, -166689119, -70229278, -60056229, -69151756, +52212618, -9135853, -64606639, +24540856, +25584363, -69501722 },
};

static int factor_smooth(ulong * c, ulong n)
{
    slong i;

    for (i = 0; i < NUM_SMALL_PRIMES; i++)
        c[i] = 0;

    /* Hardcoded so that the compiler can remove divisions */
    while (n != 1 && n % 2 == 0) { n /= 2; c[0]++; }
    while (n != 1 && n % 3 == 0) { n /= 3; c[1]++; }
    while (n != 1 && n % 5 == 0) { n /= 5; c[2]++; }
    while (n != 1 && n % 7 == 0) { n /= 7; c[3]++; }
    while (n != 1 && n % 11 == 0) { n /= 11; c[4]++; }
    while (n != 1 && n % 13 == 0) { n /= 13; c[5]++; }
    while (n != 1 && n % 17 == 0) { n /= 17; c[6]++; }
    while (n != 1 && n % 19 == 0) { n /= 19; c[7]++; }
    while (n != 1 && n % 23 == 0) { n /= 23; c[8]++; }
    while (n != 1 && n % 29 == 0) { n /= 29; c[9]++; }
    while (n != 1 && n % 31 == 0) { n /= 31; c[10]++; }
    while (n != 1 && n % 37 == 0) { n /= 37; c[11]++; }
    while (n != 1 && n % 41 == 0) { n /= 41; c[12]++; }

    return n == 1;
}

FLINT_TLS_PREFIX arb_struct _arb_log_p_cache[NUM_SMALL_PRIMES];
FLINT_TLS_PREFIX slong _arb_log_p_cache_prec = 0;

void _arb_log_p_cleanup(void)
{
    slong i;
    for (i = 0; i < NUM_SMALL_PRIMES; i++)
        arb_clear(_arb_log_p_cache + i);
    _arb_log_p_cache_prec = 0;
}

static void atanh_bs(arb_t s, ulong p, ulong q, slong prec);

typedef struct
{
    arb_ptr xs;
    slong prec;
}
log_p_precomp_work;

static void
log_p_precomp_worker(slong i, log_p_precomp_work * work)
{
    atanh_bs(work->xs + i, 1, small_primes_atanh_x[i], work->prec);
}

void _arb_log_p_ensure_cached(slong prec)
{
    slong i, wp;

    if (_arb_log_p_cache_prec < prec)
    {
        arb_ptr xs;
        log_p_precomp_work work;

        if (_arb_log_p_cache_prec == 0)
        {
            for (i = 0; i < NUM_SMALL_PRIMES; i++)
                arb_init(_arb_log_p_cache + i);

            flint_register_cleanup_function(_arb_log_p_cleanup);
        }

        wp = prec + 32;

        if (wp <= ARB_LOG_TAB2_PREC - 16)
        {
            for (i = 0; i < NUM_SMALL_PRIMES; i++)
            {
                slong exp, exp_fix;
                mp_size_t n;
                arb_ptr res = _arb_log_p_cache + i;

                n = ARB_LOG_TAB2_PREC / FLINT_BITS;

                /* exponent of log(prime(i+1)) */
                exp = (i >= 1) + (i >= 4) + (i >= 16) + (i >= 429);

                /* just reading the table is known to give the correct rounding */
                _arf_set_round_mpn(arb_midref(res), &exp_fix, arb_log_p_tab[i], n, 0, wp, ARF_RND_NEAR);
                exp += exp_fix;
                _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp);

                /* 1/2 ulp error */
                _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), exp - wp);
                MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
            }
        }
        else
        {
            prec = FLINT_MAX(prec, _arb_log_p_cache_prec * 1.25);
            wp = prec + 32;

            xs = _arb_vec_init(NUM_SMALL_PRIMES);
            work.xs = xs;
            work.prec = wp;
            flint_parallel_do((do_func_t) log_p_precomp_worker, &work, NUM_SMALL_PRIMES, -1, FLINT_PARALLEL_STRIDED);

            for (i = 0; i < NUM_SMALL_PRIMES; i++)
            {
                arb_dot_si(_arb_log_p_cache + i, NULL, 1, xs, 1, small_primes_atanh_c[i], 1, NUM_SMALL_PRIMES, wp);
                arb_mul_2exp_si(_arb_log_p_cache + i, _arb_log_p_cache + i, 1);
            }

            _arb_vec_clear(xs, NUM_SMALL_PRIMES);
        }

        _arb_log_p_cache_prec = prec;
    }
}

static void
bsplit(fmpz_t p1, fmpz_t q1, fmpz_t r1,
        const fmpz_t p, const fmpz_t q, slong a, slong b, int cont)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            fmpz_set(p1, p);
            fmpz_mul_ui(q1, q, 2*a+1);
            fmpz_mul_ui(r1, p, 2*a+1);
        }
        else
        {
            fmpz_mul(p1, p, p);
            fmpz_mul(q1, q, q);
            fmpz_mul_ui(q1, q1, 2*a+1);
            fmpz_mul_ui(r1, p1, 2*a+1);
        }
    }
    else
    {
        fmpz_t p2, q2, r2;
        slong m;

        m = (a + b) / 2;

        fmpz_init(p2);
        fmpz_init(q2);
        fmpz_init(r2);

        bsplit(p1, q1, r1, p, q, a, m, 1);
        bsplit(p2, q2, r2, p, q, m, b, cont);

        fmpz_mul(p1, p1, q2);
        fmpz_addmul(p1, r1, p2);
        fmpz_mul(q1, q1, q2);
        if (cont)
            fmpz_mul(r1, r1, r2);

        fmpz_clear(p2);
        fmpz_clear(q2);
        fmpz_clear(r2);

    }
}

#define LOG2 0.69314718055994530942

/* Assumption: p/q <= 2 */
static void
atanh_bs(arb_t s, ulong p, ulong q, slong prec)
{
    fmpz_t pp, qq, P, Q, R;
    double logqp;
    slong N;

    if (p == 0)
    {
        arb_zero(s);
        return;
    }

    fmpz_init(pp);
    fmpz_init(qq);
    fmpz_init(P);
    fmpz_init(Q);
    fmpz_init(R);

    /* If N >= 1 and p/q <= 1/2, the error is bounded by (p/q)^(2N+1).
    For error <= 2^-prec, it is sufficient to pick
    N >= (1/2) * (prec * log(2) / log(q/p) - 1). */
    logqp = mag_d_log_lower_bound(q / (double) p) * (1.0 - 1e-12);
    N = ceil((prec * (0.5 * LOG2) / logqp) * (1.0 + 1e-12));
    N = FLINT_MAX(N, 1);

    fmpz_set_ui(pp, p);
    fmpz_set_ui(qq, q);

    bsplit(P, Q, R, pp, qq, 0, N, 0);

    arb_fmpz_div_fmpz(s, P, Q, prec);
    arb_add_error_2exp_si(s, -prec);

    fmpz_clear(pp);
    fmpz_clear(qq);
    fmpz_clear(P);
    fmpz_clear(Q);
    fmpz_clear(R);
}

static int n_width(ulong k)
{
    int a, b;
    count_leading_zeros(a, k);
    count_trailing_zeros(b, k);
    return FLINT_BITS - a - b;
}

void
arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, slong prec)
{
    if (prev < 2 || prec < 600 ||
        (prec < ARB_LOG_TAB2_PREC - 64 && n_width(k) <= ARB_LOG_TAB21_BITS + 1)
        || k < prev || (k + prev) < prev ||
        (k - prev) >= 0.25 * (k + prev))
    {
        arf_t t;
        arf_init_set_ui(t, k);
        arb_log_arf(s, t, prec);
        /* no need to clear t */
    }
    else
    {
        arb_t t;
        ulong p, q;

        arb_init(t);

        p = k - prev;
        q = k + prev;

        if ((p % 2 == 0) && (q % 2 == 0))
        {
            p >>= 1;
            q >>= 1;
        }

        atanh_bs(t, p, q, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_add(s, log_prev, t, prec);

        arb_clear(t);
    }
}

int
_arb_log_ui_smooth(arb_t res, ulong n, slong prec)
{
    ulong c[NUM_SMALL_PRIMES];

    if (factor_smooth(c, n))
    {
        _arb_log_p_ensure_cached(prec);
        arb_dot_ui(res, NULL, 0, _arb_log_p_cache, 1, c, 1, NUM_SMALL_PRIMES, prec);
        return 1;
    }
    else
    {
        return 0;
    }
}

void
arb_log_ui(arb_t z, ulong x, slong prec)
{
    if (x == 2)
    {
        arb_const_log2(z, prec);
    }
    else if (x == 10)
    {
        arb_const_log10(z, prec);
    }
    else
    {
        arf_t t;
        arf_init(t);
        arf_set_ui(t, x);
        arb_log_arf(z, t, prec);
        arf_clear(t);
    }
}

void
arb_log_fmpz(arb_t z, const fmpz_t x, slong prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_log_arf(z, t, prec);
    arf_clear(t);
}
