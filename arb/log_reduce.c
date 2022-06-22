/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define TERMINATOR -32768

/* Precomputed data for ARB_LOG_PRIME_CACHE_NUM primes */

static const short small_primes[ARB_LOG_PRIME_CACHE_NUM] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41 };

/* log2 of primes */
static const float log_weights[ARB_LOG_PRIME_CACHE_NUM] = {
    0.0,   /* log(2) is free */
    1.5849625007210761,
    2.3219280948874257,
    2.807354922057584,
    3.4594316186371543,
    3.700439718141297,
    4.087462841250272,
    4.247927513443756,
    4.523561956057165,
    4.857980995127946,
    4.954196310386578,
    5.209453365629088,
    5.357552004617901,
};

/* Output of _arb_log_precompute_reductions with n = 13, C = 8 */

static const short log_rel_d[] = {
    1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0,
    1, -1, -1, 0, 1, 1, 0, -1, 0, 0, 0, 0, 0,
    1, 1, 0, -1, 0, 1, -1, -1, 0, 1, 0, 0, 0,
    1, 0, 0, 0, 1, 0, -1, 1, -1, 1, -1, 0, 0,
    0, 1, 1, -2, -1, -1, 0, 1, 1, -1, 1, 0, 0,
    0, 1, 0, -1, -1, 0, 2, -1, 0, -1, -1, 1, 1,
    1, -1, 0, 1, 1, 2, -1, 0, -2, 1, -1, -1, 1,
    1, 0, 4, -1, -2, 0, 0, 2, 0, -2, -2, 1, 1,
    0, -2, 0, 0, -2, 0, 0, 2, -4, 4, -1, 1, 0,
    1, 1, 4, 1, -1, 1, -2, -3, 0, -4, 3, 1, 1,
    3, 5, 1, 1, 0, -2, -2, 1, 1, 0, -4, 1, 1,
    0, -2, -1, 0, 2, 4, 4, 0, 3, 1, -6, -1, -3,
    3, 2, -1, -6, 2, 3, -2, -2, 3, 1, 5, -4, -2,
    0, -5, -1, -1, -9, 7, 2, -9, 3, 7, -3, 3, 0,
    1, -1, -7, -2, 5, 5, -6, 2, 0, -10, 5, 2, 3,
    3, -2, -7, -9, 6, 6, 3, 9, 1, 8, -15, -4, 0,
    10, -7, -4, -17, 4, -3, 7, -10, 9, 1, 7, 2, -4,
    2, -3, 2, -2, 15, -16, -4, 7, -11, 15, 0, -9, 4,
    1, -13, 6, 2, 2, 0, 17, -2, -11, -21, 17, 2, -2,
    6, -9, 0, 9, 9, -2, -4, -22, 4, -7, 0, 5, 11,
    2, 27, -5, 3, 15, 18, -16, 9, 0, -15, -9, -4, 1,
    13, 9, -30, 7, 9, -23, -3, 16, -21, 17, -5, 5, 6,
    4, -5, 8, -8, 6, -25, -38, -16, 24, 13, -10, 10, 24,
    36, -49, 11, 12, -24, -25, 7, 7, 7, 44, -44, 19, -5,
    101, -32, -31, 6, 36, -11, 30, 20, 8, -11, -47, -7, 0,
    48, -31, 21, -27, 34, -23, -29, 41, -50, -65, 33, 20, 40,
    58, -85, 42, -67, 51, 37, 61, 4, -65, -40, -9, 18, 5,
    68, 31, 2, 50, -12, 27, -42, 147, -38, -81, -7, -24, -9,
    26, 20, -35, 16, -1, 75, -13, 2, -128, -100, 130, 46, -13,
    0, 78, 91, -120, 149, -170, 156, 22, 5, 54, -80, -105, 8,
    226, -31, -128, -57, -35, -20, 46, -103, -5, 176, -118, 70, 21,
    137, -26, 127, 45, -14, -73, -66, -166, 71, 76, 122, -154, 53,
    146, -222, 350, -77, -31, 20, -49, -162, -179, 116, 43, 25, 81,
    286, 326, -203, 64, 71, 158, 60, -159, -69, 53, -194, 17, 3,
    254, -11, -64, 83, -168, -28, 302, -225, 459, -403, -77, 37, 29,
    169, -89, 290, -118, -605, 1100, -34, 112, 384, -141, -468, -266, -6,
    145, 375, 543, -140, -101, -298, 451, -172, -839, -201, 670, 92, -55,
    29, -159, -25, -881, 467, 368, 356, 476, -238, -301, -15, 249, -445,
    200, -51, -52, -515, 1004, -73, 319, -303, -239, -352, 674, -176, -262,
    201, 211, 1009, 127, -518, -180, 263, 108, -941, 742, -219, -722, 595,
    355, 1305, -493, -1371, 230, -720, 463, 99, 1334, 238, -612, -129, -254,
    615, 881, -82, -16, -1392, 229, 1458, 887, -961, 579, -1042, -125, -35,
    467, 771, 929, -1054, 222, -1676, 1466, -148, 642, 287, 176, -344, -783,
    556, -1642, 158, -934, 1001, -19, -1652, -329, 280, 1866, -1196, 1141, -241,
    820, -1599, 612, -2431, 2075, 187, 265, -1984, -1622, 227, 352, 1459, 650,
    2179, -1220, -281, 1277, 1083, -796, 95, -1410, -1120, 405, 2450, -3100, 1630,
    747, -2433, 64, -671, -1482, -1592, 921, 1140, 1339, -723, 3483, -1933, -462,
    1432, 1213, -345, 1948, 2565, 796, -826, -2550, -2459, 1128, -1033, -1167, 2092,
    1221, 2993, -4115, 5545, -506, -1116, 1103, -2881, -2076, -955, 2418, -675, 1344,
    1552, -1319, 1811, 403, -1630, -2912, 3300, 4551, -3903, -56, -2748, 1501, 470,
    3704, -2349, 1775, 3414, 1348, -2094, -3296, 3503, -8959, -260, 6149, 742, -849,
    5687, -1783, -4739, -4752, 2449, 2663, -640, -3687, -3290, -835, 9424, -3052, 1789,
    12240, 4587, -3750, 1958, -2949, 4395, -5453, 9129, -538, 5578, -8137, -3979, -462,
    4747, 6252, -10020, 14096, -3541, 2170, 15553, -8762, -12673, 1007, 3368, -159, -3083,
    957, 1485, -4231, 23117, -5115, -5577, 3009, -5806, -726, -8789, 1836, 17445, -11513,
    5052, -1053, -26095, -2312, -12870, 13692, 994, 5315, 12997, 9530, -7231, -3893, -3373,
    1708, 10103, -13905, 43, -1923, 3702, -4720, -19980, -10488, 6503, 29007, 3602, -6543,
    TERMINATOR,
};

static const double log_rel_epsilon[] = {
    0.18232155679395462,
    0.0266682470821613254,
    0.00350263055120206366,
    0.000442184398095918677,
    8.24980407183236868e-05,
    9.84232593854523219e-06,
    1.52063732951784046e-06,
    3.23154886181859914e-08,
    4.38246503661907078e-09,
    -2.11704006951314627e-10,
    -7.07426669244763431e-11,
    -2.21701523380274833e-11,
    3.33036400488437696e-12,
    2.54265619979145741e-13,
    -7.00024276392618185e-14,
    -9.51714644892047648e-15,
    6.80693234263194677e-16,
    1.38172686218446643e-16,
    -8.19305834374777461e-18,
    3.5176906074023809e-18,
    4.67118192241310654e-19,
    8.17891052982097575e-20,
    -1.19179263459574717e-20,
    -8.51393042738715997e-23,
    3.76392698342906217e-23,
    5.28766891158120376e-24,
    5.2060634410123273e-26,
    -1.63504294042364847e-26,
    1.90067819717850597e-27,
    -3.33141349620465382e-29,
    2.07178082822512832e-30,
    -4.74225218175538935e-31,
    -1.42274363523616913e-31,
    -7.47394389255809103e-33,
    -1.73829895124503976e-33,
    8.26424285905383777e-36,
    -4.71875282737164478e-37,
    4.63703613472910643e-37,
    2.80042458117148945e-38,
    -3.73163451438722246e-39,
    -1.74597092800142613e-40,
    -1.88295877456676773e-40,
    -2.26444799388392757e-41,
    5.60832045421886576e-43,
    -1.00422711325320322e-43,
    2.92831247865346646e-44,
    -1.13961060638254883e-45,
    -1.32337722847336868e-45,
    1.83766622090819767e-46,
    5.95732915783458014e-48,
    -2.24487689239035057e-48,
    -1.56518141370651064e-49,
    -2.78642033949683251e-51,
    5.88762413433732367e-52,
    3.40338106058818415e-53,
    -4.48648864933766729e-54,
    1.23840947437945835e-54,
    -8.98741915925760203e-56,
};

static const double log_rel_epsilon_inv[] = {
    5.48481494774707734,
    37.4977776724181737,
    285.499708114180407,
    2261.49996315130011,
    12121.4999931251623,
    101601.999999179796,
    657618.999999873224,
    30944913.5,
    228182082.833333343,
    -4723576159,
    -14135740755.5412464,
    -45105689160.5,
    300267477829.264465,
    3932895057074.63574,
    -14285218866311.6641,
    -105073511831209.594,
    1469090553077751.25,
    7237320394994939,
    -122054543986387264,
    284277417091676640,
    2.1407858152598039e+18,
    1.22265672958998426e+19,
    -8.39072143065557074e+19,
    -1.17454565612053038e+22,
    2.65679967864033022e+22,
    1.8911925400808879e+23,
    1.92083713794611083e+25,
    -6.11604732375343398e+25,
    5.26127990253408999e+26,
    -3.00172884914844711e+28,
    4.82676539128267243e+29,
    -2.10870270427044347e+30,
    -7.02867315821099744e+30,
    -1.33798167925198604e+32,
    -5.75275040742422151e+32,
    1.21003220386300296e+35,
    -2.11920402823260834e+36,
    2.15654994040373068e+36,
    3.57088709591912793e+37,
    -2.67979084271121742e+38,
    -5.72747222741376142e+39,
    -5.31079072737575268e+39,
    -4.4160872879435118e+40,
    1.78306501592245647e+42,
    -9.95790680019075144e+42,
    3.41493610155919096e+43,
    -8.77492710579701471e+44,
    -7.55642441538447421e+44,
    5.44168461400888913e+45,
    1.67860457850458665e+47,
    -4.45458725772350707e+47,
    -6.38903574526800067e+48,
    -3.58883398109482076e+50,
    1.69847798905483982e+51,
    2.93825458330304179e+52,
    -2.22891458813258854e+53,
    8.07487362369445356e+53,
    -1.11266647552533214e+55,
};

void
_arb_log_reduce_fixed(slong * rel, const short * d, const double * epsilon, const double * epsilon_inv,
    const fmpz * alpha, const float * weights,
    slong num_alpha, const fmpz_t x, slong prec, double max_weight)
{
    slong i, j, n;
    slong * new_rel;
    const short * d_row;
    double dalpha;
    double weight, dx;
    fmpz_t t;

    new_rel = flint_malloc(num_alpha * sizeof(slong));
    fmpz_init(t);

    for (i = 0; i < num_alpha; i++)
        rel[i] = 0;

    dx = fmpz_get_d(x);
    dx = ldexp(dx, -prec);

    /* Reduce by the first alpha, which is assumed to be free */
    dalpha = fmpz_get_d(alpha + 0);
    dalpha = ldexp(dalpha, -prec);
    n = floor(dx / dalpha + 0.5);
    dx -= dalpha * n;
    rel[0] = n;

    /* Recompute accurately if there is significant cancellation */
    if (FLINT_ABS(n) > 10)
    {
        fmpz_set(t, x);
        fmpz_submul_si(t, alpha + 0, rel[0]);
        dx = fmpz_get_d(t);
        dx = ldexp(dx, -prec);
    }

    for (i = 0; ; i++)
    {
        d_row = d + i * num_alpha;

        if (d_row[0] == TERMINATOR)
            break;

        for (j = 0; j < num_alpha; j++)
            new_rel[j] = d_row[j];

        n = floor(dx * epsilon_inv[i] + 0.5);

        if (n != 0)
        {
            weight = 0.0;
            for (j = 0; j < num_alpha; j++)
            {
                new_rel[j] = rel[j] + n * new_rel[j];
                if (j != 0)
                    weight += FLINT_ABS(new_rel[j]) * weights[j] * 1.442695;
            }

            if (weight > max_weight)
                break;

            for (j = 0; j < num_alpha; j++)
                rel[j] = new_rel[j];

            dx -= n * epsilon[i];
        }

        if (i % 8 == 7)
        {
            fmpz_set(t, x);
            for (j = 0; j < num_alpha; j++)
                fmpz_submul_si(t, alpha + j, rel[j]);

            dx = fmpz_get_d(t);
            dx = ldexp(dx, -prec);
        }
    }

    fmpz_clear(t);
    flint_free(new_rel);
}

static void
rel_product(fmpz_t p, fmpz_t q, const short * primes, const slong * rel, slong len)
{
    slong i;

    if (len <= 4)
    {
        fmpz_t r;
        fmpz_init(r);

        for (i = 0; i < len; i++)
        {
            fmpz_ui_pow_ui(r, primes[i], FLINT_ABS(rel[i]));

            if (rel[i] >= 0)
                fmpz_mul(p, p, r);
            else
                fmpz_mul(q, q, r);
        }

        fmpz_clear(r);
    }
    else
    {
        fmpz_t p2, q2;

        fmpz_init_set_ui(p2, 1);
        fmpz_init_set_ui(q2, 1);

        rel_product(p, q, primes, rel, len / 2);
        rel_product(p2, q2, primes + len / 2, rel + len / 2, len - len / 2);

        fmpz_mul(p, p, p2);
        fmpz_mul(q, q, q2);

        fmpz_clear(p2);
        fmpz_clear(q2);
    }
}

/* todo: error propagation */
void
_arb_exp_arf_precomp(arb_t res, const arf_t x, slong prec, int minus_one,
    slong num_logs, arb_srcptr logs, const short * primes,
    const float * weights,
    const short * log_rel_d,
    const double * epsilon, const double * epsilon_inv, double max_weight)
{
    arb_t t;
    fmpz_t p, r, q;
    slong wp;
    slong * rel;
    slong i;
    fmpz * alpha;
    slong mag;
    mag_t err, err2;

    mag = arf_abs_bound_lt_2exp_si(x);

    arb_init(t);

    rel = flint_malloc(num_logs * sizeof(slong));

    alpha = _fmpz_vec_init(num_logs);
    fmpz_init(r);

    if (prec <= 10000)
        wp = 256;
    else if (prec <= 100000)
        wp = 512;
    else
        wp = 768;

    for (i = 0; i < num_logs; i++)
        arf_get_fmpz_fixed_si(alpha + i, arb_midref(logs + i), -wp);

    arf_get_fmpz_fixed_si(r, x, -wp);
    _arb_log_reduce_fixed(rel, log_rel_d, epsilon, epsilon_inv,
        alpha, weights, num_logs, r, wp, max_weight);

    fmpz_clear(r);
    _fmpz_vec_clear(alpha, num_logs);

    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);
    if (minus_one && mag < 0)
        wp += (-mag);
    else if (mag > 0)
        wp += mag;

    arb_set_arf(t, x);
    arb_dot_si(t, t, 1, logs, 1, rel, 1, num_logs, wp);
    arb_exp_arf_generic(res, arb_midref(t), wp, 0);

    /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) */
    mag_init(err);
    mag_init(err2);
    arb_get_mag(err, res);
    mag_expm1(err2, arb_radref(t));
    mag_mul(arb_radref(res), err, err2);
    mag_clear(err);
    mag_clear(err2);

    fmpz_init(p);
    fmpz_init(q);

    fmpz_one(p);
    fmpz_one(q);
    rel_product(p, q, primes + 1, rel + 1, num_logs - 1);

    arb_mul_fmpz(res, res, p, wp);
    arb_div_fmpz(res, res, q, wp);
    arb_mul_2exp_si(res, res, rel[0]);

    if (minus_one)
        arb_sub_ui(res, res, 1, prec);
    else
        arb_set_round(res, res, prec);

    flint_free(rel);

    fmpz_clear(p);
    fmpz_clear(q);
    arb_clear(t);
}

void arb_exp_arf_huge(arb_t z, const arf_t x, slong mag, slong prec, int minus_one);

void
arb_exp_arf_log_reduction(arb_t res, const arf_t x, slong prec, int minus_one)
{
    slong wp;
    slong mag;

    mag = arf_abs_bound_lt_2exp_si(x);

    if (mag < -prec / 16 || mag < -768 || arf_bits(x) < prec / 128)
    {
        arb_exp_arf_generic(res, x, prec, minus_one);
        return;
    }

    /* multiprecision log(2) reduction not implemented here */
    if ((FLINT_BITS == 32 && mag > 20) || (FLINT_BITS == 64 && mag > 40))
    {
        arb_exp_arf_huge(res, x, mag, prec, minus_one);
        return;
    }

    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);
    wp += FLINT_BITS;
    if (minus_one && mag < 0)
        wp += (-mag);
    else if (mag > 0)
        wp += mag;

    _arb_log_p_ensure_cached(wp);

    _arb_exp_arf_precomp(res, x, prec, minus_one,
        ARB_LOG_PRIME_CACHE_NUM,
        _arb_log_p_cache,
        small_primes, log_weights,
        log_rel_d,
        log_rel_epsilon, log_rel_epsilon_inv, prec);
}
