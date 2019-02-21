/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_hypgeom.h"
#include "acb_calc.h"
#include "acb_dft.h"

void platt_lemma_32(arb_t out, const arb_t h, const arb_t t0,
        const arb_t x, slong prec)
{
    arb_t pi, one_fourth;
    arb_t x1, x2;

    arb_init(pi);
    arb_init(one_fourth);
    arb_init(x1);
    arb_init(x2);

    arb_const_pi(pi, prec);
    arb_set_d(one_fourth, 0.25);

    arb_set_d(x1, 1.25);
    arb_pow(x1, pi, x1, prec);
    arb_mul_2exp_si(x1, x1, 1);

    arb_sqr(x2, t0, prec);
    arb_sub(x2, x2, one_fourth, prec);
    arb_div(x2, x2, h, prec);
    arb_div(x2, x2, h, prec);
    arb_mul_2exp_si(x2, x2, -1);

    arb_mul(out, pi, x, prec);
    arb_add(out, out, x2, prec);
    arb_neg(out, out);
    arb_exp(out, out, prec);
    arb_mul(out, out, x1, prec);

    arb_clear(pi);
    arb_clear(one_fourth);
    arb_clear(x1);
    arb_clear(x2);
}

void platt_lemma_A5(arb_t out,
        const arb_t B, const arb_t h, slong k, slong prec)
{
    arb_t a, b;
    arb_t x1, x2, x3, x4, x5, x6;

    arb_init(a);
    arb_init(b);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);
    arb_init(x6);

    arb_const_pi(x1, prec);
    arb_mul(x1, x1, B, prec);
    arb_pow_ui(x1, x1, (ulong) k, prec);
    arb_mul_2exp_si(x1, x1, 3);

    arb_div(a, B, h, prec);
    arb_sqr(a, a, prec);
    arb_mul_2exp_si(a, a, -3);

    arb_neg(x2, a);
    arb_exp(x2, x2, prec);

    arb_set_si(x3, 3*k - 1);
    arb_mul_2exp_si(x3, x3, -1);

    arb_set_ui(x4, 2);
    arb_pow(x4, x4, x3, prec);

    arb_set_si(b, k+1);

    arb_div(x5, h, B, prec);
    arb_pow_ui(x5, x5, (ulong) (k+1), prec);

    arb_mul_2exp_si(x6, b, -1);

    arb_hypgeom_gamma_upper(x6, x6, a, 0, prec);

    arb_mul(out, x4, x5, prec);
    arb_mul(out, out, x6, prec);
    arb_add(out, out, x2, prec);
    arb_mul(out, out, x1, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
    arb_clear(x6);
}

void platt_lemma_A3_X(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k, slong prec)
{
    slong l;
    arb_t total, summand;
    slong upper_limit;
    arb_t two, half;
    arb_t x1, x2, x3, x4, x5, x6;
    arb_t a;

    arb_init(total);
    arb_init(summand);
    arb_init(two);
    arb_init(half);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);
    arb_init(x6);
    arb_init(a);

    arb_set_ui(two, 2);
    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    upper_limit = (sigma-1)/2;

    arb_add_si(a, half, sigma, prec);
    arb_div(a, a, h, prec);
    arb_sqr(a, a, prec);
    arb_mul_2exp_si(a, a, -1);

    for (l=0; l<=upper_limit; l++)
    {
        arb_bin_uiui(x1, (ulong) upper_limit, (ulong) l, prec);

        arb_pow_ui(x2, t0, (ulong) (upper_limit - l), prec);

        arb_pow_ui(x3, h, (ulong) (k + l + 1), prec);

        arb_set_si(x4, k+l-1);
        arb_mul_2exp_si(x4, x4, -1);
        arb_pow(x4, two, x4, prec);

        arb_set_si(x5, k + l + 1);
        arb_mul_2exp_si(x5, x5, -1);

        arb_hypgeom_gamma_upper(x6, x5, a, 0, prec);

        arb_mul(summand, x1, x2, prec);
        arb_mul(summand, summand, x3, prec);
        arb_mul(summand, summand, x4, prec);
        arb_mul(summand, summand, x6, prec);

        arb_add(total, total, summand, prec);
    }

    arb_set(out, total);

    arb_clear(total);
    arb_clear(summand);
    arb_clear(two);
    arb_clear(half);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
    arb_clear(x6);
    arb_clear(a);
}

void platt_lemma_A3(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k, slong prec)
{
    arb_t two, half, pi, X;
    arb_t x1, x2, x3, x4, x5, x6, x7;
    arb_t z1, z2;

    arb_init(two);
    arb_init(half);
    arb_init(pi);
    arb_init(X);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);
    arb_init(x6);
    arb_init(x7);
    arb_init(z1);
    arb_init(z2);

    arb_set_ui(two, 2);
    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_const_pi(pi, prec);

    arb_set_si(x1, 6*k + 5 - sigma);
    arb_mul_2exp_si(x1, x1, -2);
    arb_pow(x1, two, x1, prec);

    arb_pow_ui(x2, pi, (ulong) k, prec);

    arb_sqrt_ui(x3, 2, prec);
    arb_mul_2exp_si(x3, x3, 1);
    arb_add_ui(x3, x3, 1, prec);
    arb_div_ui(x3, x3, 6, prec);
    arb_div(x3, x3, t0, prec);
    arb_exp(x3, x3, prec);

    arb_add_si(x4, half, sigma, prec);
    arb_pow_ui(x4, x4, (ulong) k, prec);

    arb_add_si(x5, half, sigma, prec);
    arb_add(x5, x5, t0, prec);
    arb_pow_ui(x5, x5, (ulong) ((sigma - 1)/2), prec);

    arb_mul(z1, x1, x2, prec);
    arb_mul(z1, z1, x3, prec);
    arb_mul(z1, z1, x4, prec);
    arb_mul(z1, z1, x5, prec);
    arb_mul(z1, z1, h, prec);

    arb_set_si(x6, 6*k + 7 - sigma);
    arb_mul_2exp_si(x6, x6, -2);
    arb_pow(x6, two, x6, prec);

    arb_set_si(x7, 2*k-1);
    arb_mul_2exp_si(x7, x7, -1);
    arb_pow(x7, pi, x7, prec);

    platt_lemma_A3_X(X, sigma, t0, h, k, prec);

    arb_mul(z2, x6, x7, prec);
    arb_mul(z2, z2, x3, prec);
    arb_mul(z2, z2, X, prec);

    arb_add(out, z1, z2, prec);

    arb_clear(two);
    arb_clear(half);
    arb_clear(pi);
    arb_clear(X);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
    arb_clear(x6);
    arb_clear(x7);
    arb_clear(z1);
    arb_clear(z2);
}

void platt_lemma_A7_S(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k,
        const arb_t A, slong prec)
{
    slong l;
    arb_t total, summand;
    arb_t pi, half;
    arb_t a;
    arb_t l_factorial, kd2, t02;
    arb_t x1, x2, x3, x4, x5;

    arb_init(total);
    arb_init(summand);
    arb_init(pi);
    arb_init(half);
    arb_init(a);
    arb_init(l_factorial);
    arb_init(kd2);
    arb_init(t02);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_const_pi(pi, prec);
    arb_one(l_factorial);
    arb_set_si(kd2, k);
    arb_mul_2exp_si(kd2, kd2, -1);
    arb_sqr(t02, t0, prec);

    for (l=0; l<=(sigma-1)/2; l++)
    {
        if (l > 1)
        {
            arb_mul_si(l_factorial, l_factorial, l, prec);
        }

        arb_mul_si(a, pi, 4*l+1, prec);
        arb_mul(a, a, A, prec);

        arb_inv(x1, a, prec);
        arb_add_ui(x1, x1, 1, prec);

        arb_add_ui(x2, half, (ulong) (2*l), prec);
        arb_sqr(x2, x2, prec);
        arb_add(x2, x2, t02, prec);
        arb_pow(x2, x2, kd2, prec);
        arb_div(x2, x2, l_factorial, prec);

        arb_set_si(x3, 4*l + 1);
        arb_div(x3, x3, h, prec);
        arb_sqr(x3, x3, prec);
        arb_mul_2exp_si(x3, x3, -3);

        arb_mul_2exp_si(x4, a, -1);

        arb_sub(x5, x3, x4, prec);
        arb_exp(x5, x5, prec);

        arb_mul(summand, x1, x2, prec);
        arb_mul(summand, summand, x5, prec);

        arb_add(total, total, summand, prec);
    }

    arb_set(out, total);

    arb_clear(total);
    arb_clear(summand);
    arb_clear(pi);
    arb_clear(half);
    arb_clear(a);
    arb_clear(l_factorial);
    arb_clear(kd2);
    arb_clear(t02);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

void platt_lemma_A7(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k,
        const arb_t A, slong prec)
{
    arb_t S, C;
    arb_t pi, a;
    arb_t x1, x2;
    arb_t y1, y2, y3, y4;
    arb_t z1, z2;

    arb_init(S);
    arb_init(C);
    arb_init(pi);
    arb_init(a);
    arb_init(x1);
    arb_init(x2);
    arb_init(y1);
    arb_init(y2);
    arb_init(y3);
    arb_init(y4);
    arb_init(z1);
    arb_init(z2);

    arb_const_pi(pi, prec);

    arb_pow_ui(x1, pi, (ulong) (k+1), prec);
    arb_mul_2exp_si(x1, x1, k+3);

    arb_div(x2, t0, h, prec);
    arb_sqr(x2, x2, prec);
    arb_mul_2exp_si(x2, x2, -1);
    arb_neg(x2, x2);
    arb_exp(x2, x2, prec);

    platt_lemma_A7_S(S, sigma, t0, h, k, A, prec);

    arb_mul(z1, x1, x2, prec);
    arb_mul(z1, z1, S, prec);

    arb_mul_si(a, pi, 2*sigma-1, prec);
    arb_mul(a, a, A, prec);

    arb_inv(y1, a, prec);
    arb_add_ui(y1, y1, 1, prec);

    arb_set_si(y2, 2*sigma + 1);
    arb_div(y2, y2, h, prec);
    arb_sqr(y2, y2, prec);
    arb_mul_2exp_si(y2, y2, -3);

    arb_mul_2exp_si(y3, a, -1);

    arb_sub(y4, y2, y3, prec);
    arb_exp(y4, y4, prec);

    platt_lemma_A3(C, sigma, t0, h, k, prec);

    arb_mul(z2, y1, y4, prec);
    arb_mul(z2, z2, C, prec);
    arb_mul_2exp_si(z2, z2, 1);

    arb_add(out, z1, z2, prec);

    arb_clear(S);
    arb_clear(C);
    arb_clear(pi);
    arb_clear(a);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(y1);
    arb_clear(y2);
    arb_clear(y3);
    arb_clear(y4);
    arb_clear(z1);
    arb_clear(z2);
}

void platt_lemma_A9_a(arb_t out, slong sigma,
        const arb_t t0, const arb_t h,
        const arb_t A, slong prec)
{
    arb_t a, pi, C;
    arb_t y1, y2, y3, y4;
    arb_t z1, z2, z3;
    
    arb_init(a);
    arb_init(pi);
    arb_init(C);
    arb_init(y1);
    arb_init(y2);
    arb_init(y3);
    arb_init(y4);
    arb_init(z1);
    arb_init(z2);
    arb_init(z3);

    arb_const_pi(pi, prec);

    arb_mul_si(a, pi, 2*sigma-1, prec);
    arb_mul(a, a, A, prec);

    arb_inv(y1, a, prec);
    arb_add_ui(y1, y1, 1, prec);

    arb_set_si(y2, 2*sigma - 1);
    arb_div(y2, y2, h, prec);
    arb_sqr(y2, y2, prec);
    arb_mul_2exp_si(y2, y2, -3);

    arb_mul_2exp_si(y3, a, -1);

    arb_sub(y4, y2, y3, prec);
    arb_exp(y4, y4, prec);

    platt_lemma_A3(C, sigma, t0, h, 0, prec);

    arb_zeta_ui(z1, (ulong) sigma, prec);
    arb_mul_2exp_si(z1, z1, 1);

    arb_set_si(z2, 1-2*sigma);
    arb_mul_2exp_si(z2, z2, -2);
    arb_pow(z2, pi, z2, prec);

    arb_sub(z3, y2, y3, prec);
    arb_exp(z3, z3, prec);

    arb_mul(out, z1, z2, prec);
    arb_mul(out, out, z3, prec);
    arb_mul(out, out, C, prec);
    arb_mul(out, out, y1, prec);

    arb_clear(a);
    arb_clear(pi);
    arb_clear(C);
    arb_clear(y1);
    arb_clear(y2);
    arb_clear(y3);
    arb_clear(y4);
    arb_clear(z1);
    arb_clear(z2);
    arb_clear(z3);
}

void platt_lemma_A9_b(arb_t out,
        const arb_t t0, const arb_t h, const arb_t A, slong prec)
{
    arb_t pi;
    arb_t x1, x2, x3, x4, x5;

    arb_init(pi);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_const_pi(pi, prec);

    arb_set_ui(x1, 5);
    arb_mul_2exp_si(x1, x1, -2);
    arb_pow(x1, pi, x1, prec);
    arb_mul_2exp_si(x1, x1, 2);

    arb_sqr(x2, t0, prec);
    arb_mul_2exp_si(x2, x2, 2);
    arb_sub_ui(x2, x2, 1, prec);
    arb_neg(x2, x2);
    arb_div(x2, x2, h, prec);
    arb_div(x2, x2, h, prec);
    arb_mul_2exp_si(x2, x2, -3);

    arb_mul(x3, A, pi, prec);
    arb_mul_2exp_si(x3, x3, -1);

    arb_mul(x4, A, pi, prec);
    arb_inv(x4, x4, prec);
    arb_add_ui(x4, x4, 1, prec);

    arb_sub(x5, x2, x3, prec);
    arb_exp(x5, x5, prec);

    arb_mul(out, x1, x4, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(pi);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}


void platt_lemma_A9(arb_t out, slong sigma,
        const arb_t t0, const arb_t h,
        const arb_t A, slong prec)
{
    arb_t a, b;

    arb_init(a);
    arb_init(b);

    platt_lemma_A9_a(a, sigma, t0, h, A, prec);
    platt_lemma_A9_b(b, t0, h, A, prec);

    arb_add(out, a, b, prec);

    arb_clear(a);
    arb_clear(b);
}

void platt_lemma_A11_X(arb_t out,
        const arb_t t0, const arb_t h, const arb_t B,
        const arb_t beta, slong prec)
{
    arb_t x1, x2;

    arb_init(x1);
    arb_init(x2);

    arb_mul_2exp_si(x1, B, -1);
    arb_add(x1, x1, t0, prec);
    arb_pow(x1, x1, beta, prec);

    arb_div(x2, B, h, prec);
    arb_sqr(x2, x2, prec);
    arb_mul_2exp_si(x2, x2, -3);
    arb_neg(x2, x2);
    arb_exp(x2, x2, prec);

    arb_mul(out, x1, x2, prec);

    arb_clear(x1);
    arb_clear(x2);
}

void platt_lemma_A11_Y(arb_t out,
        const arb_t t0, const arb_t h, const arb_t B,
        const arb_t beta, slong prec)
{
    arb_t x1, x2, x3, x4, x5;

    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_rsqrt_ui(x1, 2, prec);

    arb_pow(x2, t0, beta, prec);

    arb_one(x3);
    arb_mul_2exp_si(x3, x3, -1);

    arb_div(x4, B, h, prec);
    arb_sqr(x4, x4, prec);
    arb_mul_2exp_si(x4, x4, -3);

    arb_hypgeom_gamma_upper(x5, x3, x4, 0, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}


void platt_lemma_A11_Z(arb_t out,
        const arb_t t0, const arb_t h,
        const arb_t beta, slong prec)
{
    arb_t two;
    arb_t x1, x2, x3, x4, x5;

    arb_init(two);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_set_ui(two, 2);

    arb_sub_ui(x1, beta, 1, prec);
    arb_mul_2exp_si(x1, x1, -1);
    arb_pow(x1, two, x1, prec);

    arb_pow(x2, h, beta, prec);

    arb_add_ui(x3, beta, 1, prec);
    arb_mul_2exp_si(x3, x3, -1);

    arb_div(x4, t0, h, prec);
    arb_sqr(x4, x4, prec);
    arb_mul_2exp_si(x4, x4, -1);

    arb_hypgeom_gamma_upper(x5, x3, x4, 0, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(two);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}


void platt_lemma_A11_beta(arb_t out, const arb_t t0, slong prec)
{
    arb_t logt0, one_sixth, beta;

    arb_init(logt0);
    arb_init(one_sixth);
    arb_init(beta);

    arb_set_ui(one_sixth, 6);
    arb_inv(one_sixth, one_sixth, prec);
    arb_log(logt0, t0, prec);

    arb_log(beta, logt0, prec);
    arb_div(beta, beta, logt0, prec);
    arb_add(beta, beta, one_sixth, prec);

    arb_set(out, beta);

    arb_clear(logt0);
    arb_clear(one_sixth);
    arb_clear(beta);
}


void platt_lemma_A11(arb_t out,
        const arb_t t0, const arb_t h, const arb_t B, slong prec)
{
    arb_t beta;
    arb_t X, Y, Z;
    arb_t x1, x2;

    arb_init(beta);
    arb_init(X);
    arb_init(Y);
    arb_init(Z);
    arb_init(x1);
    arb_init(x2);

    platt_lemma_A11_beta(beta, t0, prec);
    platt_lemma_A11_X(X, t0, h, B, beta, prec);
    platt_lemma_A11_Y(Y, t0, h, B, beta, prec);
    platt_lemma_A11_Z(Z, t0, h, beta, prec);

    arb_set_ui(x1, 2);
    arb_pow(x1, x1, beta, prec);
    arb_mul(x1, x1, h, prec);
    arb_div(x1, x1, B, prec);

    arb_add(x2, Y, Z, prec);
    arb_mul(x2, x2, x1, prec);
    arb_add(x2, x2, X, prec);
    arb_mul_ui(x2, x2, 6, prec);

    arb_set(out, x2);

    arb_clear(beta);
    arb_clear(X);
    arb_clear(Y);
    arb_clear(Z);
    arb_clear(x1);
    arb_clear(x2);
}

void platt_lemma_B1(arb_t out,
        slong sigma, const arb_t t0, const arb_t h, slong J_slong, slong prec)
{
    arb_t pi, J, C;
    arb_t x1, x2, x3;

    arb_init(pi);
    arb_init(J);
    arb_init(C);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);

    arb_const_pi(pi, prec);
    arb_set_si(J, J_slong);
    platt_lemma_A3(C, sigma, t0, h, 0, prec);

    arb_set_si(x1, 2*sigma - 1);
    arb_div(x1, x1, h, prec);
    arb_sqr(x1, x1, prec);
    arb_mul_2exp_si(x1, x1, -3);
    arb_exp(x1, x1, prec);

    arb_set_si(x2, 1 - 2*sigma);
    arb_mul_2exp_si(x2, x2, -2);
    arb_pow(x2, pi, x2, prec);

    arb_set_si(x3, 1 - sigma);
    arb_pow(x3, J, x3, prec);
    arb_div_si(x3, x3, sigma - 1, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x3, prec);
    arb_mul(out, out, C, prec);

    arb_clear(pi);
    arb_clear(J);
    arb_clear(C);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
}

void platt_lemma_B2(arb_t out,
        slong K_slong, const arb_t h, const arb_t xi, slong prec)
{
    arb_t two, half, pi, K;
    arb_t x1, x2, x3, x4, x5;

    arb_init(two);
    arb_init(half);
    arb_init(pi);
    arb_init(K);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_set_ui(two, 2);
    arb_mul_2exp_si(half, two, -2);
    arb_const_pi(pi, prec);
    arb_set_si(K, K_slong);

    arb_add_ui(x1, K, 5, prec);
    arb_mul_2exp_si(x1, x1, -1);
    arb_pow(x1, two, x1, prec);

    arb_add(x2, K, half, prec);
    arb_pow(x2, pi, x2, prec);

    arb_pow_ui(x3, h, (ulong) (K_slong + 1), prec);

    arb_pow_ui(x4, xi, (ulong) K_slong, prec);

    arb_add_ui(x5, K, 2, prec);
    arb_mul_2exp_si(x5, x5, -1);
    arb_rgamma(x5, x5, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x3, prec);
    arb_mul(out, out, x4, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(two);
    arb_clear(half);
    arb_clear(pi);
    arb_clear(K);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

void platt_g_gamma_term(acb_t out, const arb_t t0, const acb_t t, slong prec)
{
    acb_t z;
    acb_init(z);
    acb_add_arb(z, t, t0, prec);
    acb_mul_onei(z, z);
    acb_mul_2exp_si(z, z, 1);
    acb_add_ui(z, z, 1, prec);
    acb_mul_2exp_si(z, z, -2);
    acb_gamma(out, z, prec);
    acb_clear(z);
}


void platt_g_exp_term(acb_t out,
        const arb_t t0, const arb_t h, const acb_t t, slong prec)
{
    arb_t pi;
    acb_t z1, z2;
    arb_init(pi);
    acb_init(z1);
    acb_init(z2);
    arb_const_pi(pi, prec);
    acb_add_arb(z1, t, t0, prec);
    acb_mul_arb(z1, z1, pi, prec);
    acb_mul_2exp_si(z1, z1, -2);
    acb_div_arb(z2, t, h, prec);
    acb_sqr(z2, z2, prec);
    acb_mul_2exp_si(z2, z2, -1);
    acb_sub(out, z1, z2, prec);
    acb_exp(out, out, prec);
    arb_clear(pi);
    acb_clear(z1);
    acb_clear(z2);
}

static void
platt_g_base(acb_t out, const acb_t t, slong prec)
{
    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, prec);
    acb_mul_arb(out, t, pi, prec);
    acb_mul_onei(out, out);
    acb_mul_2exp_si(out, out, 1);
    acb_neg(out, out);
    arb_clear(pi);
}


void platt_g(acb_t out,
        const arb_t t0, const arb_t h, const acb_t t, slong k, slong prec)
{
    acb_t z1, z2, z3;
    acb_init(z1);
    acb_init(z2);
    acb_init(z3);
    platt_g_gamma_term(z1, t0, t, prec);
    platt_g_exp_term(z2, t0, h, t, prec);
    platt_g_base(z3, t, prec);
    acb_pow_ui(z3, z3, (ulong) k, prec);
    acb_mul(out, z1, z2, prec);
    acb_mul(out, out, z3, prec);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);
}


void platt_g_table(acb_ptr table,
        const arb_t A, const arb_t t0, const arb_t h,
        slong N, slong K, slong prec)
{
    slong i, n, k;
    acb_t t, base;
    acb_t gamma_term, exp_term, coeff;
    acb_ptr precomputed_powers;

    acb_init(t);
    acb_init(base);
    acb_init(gamma_term);
    acb_init(exp_term);
    acb_init(coeff);

    precomputed_powers = _acb_vec_init(K);

    for (i=0; i<N; i++)
    {
        n = i - N/2;

        acb_set_si(t, n);
        acb_div_arb(t, t, A, prec);

        platt_g_base(base, t, prec);
        _acb_vec_set_powers(precomputed_powers, base, K, prec);

        platt_g_gamma_term(gamma_term, t0, t, prec);
        platt_g_exp_term(exp_term, t0, h, t, prec);
        acb_mul(coeff, gamma_term, exp_term, prec);

        for (k=0; k<K; k++)
        {
            acb_mul(table + k*N + i, coeff, precomputed_powers + k, prec);
        }
    }

    acb_clear(t);
    acb_clear(base);
    acb_clear(gamma_term);
    acb_clear(exp_term);
    acb_clear(coeff);

    _acb_vec_clear(precomputed_powers, K);
}

static void
logjsqrtpi(arb_t out, slong j, slong prec)
{
    arb_const_sqrt_pi(out, prec);
    arb_mul_si(out, out, j, prec);
    arb_log(out, out, prec);
}

static double
_arb_mag_get_d_log2_approx(const arb_t x)
{
    double out;
    mag_t m;
    mag_init(m);
    arb_get_mag(m, x);
    out = mag_get_d_log2_approx(m);
    mag_clear(m);
    return out;
}

static void
_acb_add_error_arb_mag(acb_t res, const arb_t x)
{
    mag_t err;
    mag_init(err);
    arb_get_mag(err, x);
    acb_add_error_mag(res, err);
    mag_clear(err);
}

static void
_acb_vec_scalar_add_error_mag(acb_ptr res, slong len, const mag_t err)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        acb_add_error_mag(res + i, err);
    }
}

static void
_acb_vec_scalar_add_error_arb_mag(acb_ptr res, slong len, const arb_t x)
{
    mag_t err;
    mag_init(err);
    arb_get_mag(err, x);
    _acb_vec_scalar_add_error_mag(res, len, err);
    mag_clear(err);
}


slong get_smk_index(slong *out, const arb_t B, slong j, slong prec)
{
    arb_t pi, half, x;
    fmpz_t z;
    slong result = 0;

    arb_init(pi);
    arb_init(half);
    arb_init(x);
    fmpz_init(z);

    arb_const_pi(pi, prec);
    arb_set_d(half, 0.5);

    logjsqrtpi(x, j, prec);
    arb_div(x, x, pi, prec);
    arb_mul_2exp_si(x, x, -1);
    arb_mul(x, x, B, prec);
    arb_add(x, x, half, prec);
    arb_floor(x, x, prec);

    if (arb_get_unique_fmpz(z, x))
    {
        *out = fmpz_get_si(z);
        result = 1;
    }

    arb_clear(pi);
    arb_clear(half);
    arb_clear(x);
    fmpz_clear(z);

    return result;
}


slong fill_smk_index_table(slong *table, const arb_t B, slong J, slong prec)
{
    slong j;
    for (j = 1; j <= J; j++)
    {
        if (!get_smk_index(table + j - 1, B, j, prec))
        {
            return 0;
        }
    }
    return 1;
}

void hardy_z_per_unwindowed_platt_f(arb_t out, const arb_t t, slong prec)
{
    arb_t pi, half, a;
    acb_t z;

    arb_init(pi);
    arb_init(half);
    arb_init(a);
    acb_init(z);

    arb_const_pi(pi, prec);
    arb_mul(a, pi, t, prec);
    arb_mul_2exp_si(a, a, -2);
    arb_neg(a, a);

    arb_set_d(half, 0.5);
    acb_set_arb_arb(z, half, t);
    acb_mul_2exp_si(z, z, -1);
    acb_lgamma(z, z, prec);

    arb_sub(out, a, acb_realref(z), prec);
    arb_exp(out, out, prec);

    arb_clear(pi);
    arb_clear(half);
    arb_clear(a);
    acb_clear(z);
}

void hardy_z_per_platt_f(arb_t out,
        const arb_t h, const arb_t t0, slong n, slong N, slong B, slong prec)
{
    arb_t pi, half, a, b, t, x;
    acb_t z;

    arb_init(pi);
    arb_init(half);
    arb_init(a);
    arb_init(b);
    arb_init(t);
    arb_init(x);
    acb_init(z);

    prec += FLINT_MAX(0, _arb_mag_get_d_log2_approx(t0));

    arb_const_pi(pi, prec);
    arb_set_d(half, 0.5);

    arb_set_si(t, n);
    arb_mul_si(t, t, B, prec);
    arb_div_si(t, t, N, prec);

    arb_add(x, t0, t, prec);

    arb_mul(a, pi, x, prec);
    arb_mul_2exp_si(a, a, -2);
    arb_neg(a, a);

    arb_div(b, t, h, prec);
    arb_sqr(b, b, prec);
    arb_mul_2exp_si(b, b, -1);

    acb_set_arb_arb(z, half, x);
    acb_mul_2exp_si(z, z, -1);
    acb_lgamma(z, z, prec);

    arb_add(out, a, b, prec);
    arb_sub(out, out, acb_realref(z), prec);
    arb_exp(out, out, prec);

    arb_clear(pi);
    arb_clear(half);
    arb_clear(a);
    arb_clear(b);
    arb_clear(x);
    acb_clear(z);
}


slong platt_smk(acb_ptr table,
        const arb_t t0, const arb_t B, slong N, slong J, slong K, slong prec)
{
    slong result = 1;
    slong j, k, m;
    acb_ptr row;
    arb_ptr diff_powers;
    arb_t rpi, rsqrtj, um, a, base;
    acb_t z;

    arb_init(rpi);
    arb_init(rsqrtj);
    arb_init(um);
    arb_init(a);
    arb_init(base);
    acb_init(z);
    diff_powers = _arb_vec_init(K);

    arb_const_pi(rpi, prec);
    arb_inv(rpi, rpi, prec);

    for (j = 1; j <= J; j++)
    {
        logjsqrtpi(a, j, prec);
        arb_mul(a, a, rpi, prec);

        arb_rsqrt_ui(rsqrtj, (ulong) j, prec);

        acb_set_arb(z, t0);
        acb_mul_arb(z, z, a, prec);
        acb_neg(z, z);
        acb_exp_pi_i(z, z, prec);
        acb_mul_arb(z, z, rsqrtj, prec);

        if (!get_smk_index(&m, B, j, prec))
        {
            result = 0;
            goto finish;
        }
        arb_set_si(um, m);
        arb_div(um, um, B, prec);

        arb_mul_2exp_si(base, a, -1);
        arb_sub(base, base, um, prec);

        _arb_vec_set_powers(diff_powers, base, K, prec);

        for (k = 0; k < K; k++)
        {
            row = table + N*k;
            acb_addmul_arb(row + m, z, diff_powers + k, prec);
        }
    }

finish:

    arb_clear(rpi);
    arb_clear(rsqrtj);
    arb_clear(um);
    arb_clear(a);
    arb_clear(base);
    acb_clear(z);
    _arb_vec_clear(diff_powers, K);

    return result;
}

void platt_smk_precomp(acb_ptr table, const slong *smk_index_table,
        const arb_t t0, const arb_t B, slong N, slong J, slong K, slong prec)
{
    slong j, k, m;
    acb_ptr row;
    arb_ptr diff_powers;
    arb_t rpi, rsqrtj, um, a, base;
    acb_t z;

    arb_init(rpi);
    arb_init(rsqrtj);
    arb_init(um);
    arb_init(a);
    arb_init(base);
    acb_init(z);
    diff_powers = _arb_vec_init(K);

    arb_const_pi(rpi, prec);
    arb_inv(rpi, rpi, prec);

    for (j = 1; j <= J; j++)
    {
        logjsqrtpi(a, j, prec);
        arb_mul(a, a, rpi, prec);

        arb_rsqrt_ui(rsqrtj, (ulong) j, prec);

        acb_set_arb(z, t0);
        acb_mul_arb(z, z, a, prec);
        acb_neg(z, z);
        acb_exp_pi_i(z, z, prec);
        acb_mul_arb(z, z, rsqrtj, prec);

        m = smk_index_table[j - 1];
        arb_set_si(um, m);
        arb_div(um, um, B, prec);

        arb_mul_2exp_si(base, a, -1);
        arb_sub(base, base, um, prec);

        _arb_vec_set_powers(diff_powers, base, K, prec);

        for (k = 0; k < K; k++)
        {
            row = table + N*k;
            acb_addmul_arb(row + m, z, diff_powers + k, prec);
        }
    }

    arb_clear(rpi);
    arb_clear(rsqrtj);
    arb_clear(um);
    arb_clear(a);
    arb_clear(base);
    acb_clear(z);
    _arb_vec_clear(diff_powers, K);
}


void do_convolutions(acb_ptr out_table,
        const acb_ptr table, const acb_ptr S_table,
        slong N, slong K, slong prec)
{
    slong i, k;
    acb_ptr padded_table_row, padded_S_table_row, padded_out_table;
    acb_ptr fp, gp;
    acb_dft_pre_t pre;

    padded_table_row = _acb_vec_init(N*2);
    padded_S_table_row = _acb_vec_init(N*2);
    padded_out_table = _acb_vec_init(N*2);
    fp = _acb_vec_init(N*2);
    gp = _acb_vec_init(N*2);

    acb_dft_precomp_init(pre, N*2, prec);

    for (k = 0; k < K; k++)
    {
        _acb_vec_zero(padded_table_row, N*2);
        _acb_vec_zero(padded_S_table_row, N*2);
        _acb_vec_zero(padded_out_table, N*2);

        _acb_vec_set(padded_table_row, table + k*N, N);
        _acb_vec_set(padded_S_table_row, S_table + k*N, N);

        for (i = 1; i < N; i++)
        {
            acb_swap(padded_S_table_row + i, padded_S_table_row + N*2 - i);
        }

        acb_dft_precomp(fp, padded_S_table_row, pre, prec);
        acb_dft_precomp(gp, padded_table_row, pre, prec);
        _acb_vec_kronecker_mul(gp, gp, fp, N*2, prec);
        acb_dft_inverse_precomp(padded_out_table, gp, pre, prec);

        for (i = 0; i <= N/2; i++)
        {
            acb_add(out_table + i,
                    out_table + i, padded_out_table + i, prec);
        }
    }

    _acb_vec_clear(padded_table_row, N*2);
    _acb_vec_clear(padded_S_table_row, N*2);
    _acb_vec_clear(padded_out_table, N*2);
    _acb_vec_clear(fp, N*2);
    _acb_vec_clear(gp, N*2);

    acb_dft_precomp_clear(pre);
}




slong
platt_multi_evaluation(acb_ptr out,
        slong B_slong, const arb_t A, slong N, const arb_t t0, const arb_t h,
        slong J, slong K, slong sigma, slong prec)
{
    slong result = 1;
    slong i, n, k;
    acb_ptr table, S_table, out_table;
    acb_ptr row;
    arb_t B, t, x, k_factorial, err, ratio, c, xi;
    acb_t z;
    acb_dft_pre_t pre_N;

    arb_init(B);
    arb_init(t);
    arb_init(x);
    arb_init(k_factorial);
    arb_init(err);
    arb_init(ratio);
    arb_init(c);
    arb_init(xi);
    acb_init(z);
    table = _acb_vec_init(K*N);
    S_table = _acb_vec_init(K*N);
    out_table = _acb_vec_init(N);
    acb_dft_precomp_init(pre_N, N, prec);

    arb_set_si(B, B_slong);
    arb_inv(xi, B, prec);
    arb_mul_2exp_si(xi, xi, -1);

    if (!platt_smk(S_table, t0, B, N, J, K, prec))
    {
        result = 0;
        goto finish;
    }

    platt_g_table(table, A, t0, h, N, K, prec);

    for (k = 0; k < K; k++)
    {
        platt_lemma_A5(err, B, h, k, prec);
        _acb_vec_scalar_add_error_arb_mag(table + N*k, N, err);
    }

    for (k = 0; k < K; k++)
    {
        row = table + N*k;
        for (i = 0; i < N/2; i++)
        {
            acb_swap(row + i, row + i + N/2);
        }
        acb_dft_precomp(row, row, pre_N, prec);
    }
    _acb_vec_scalar_div_arb(table, table, N*K, A, prec);

    for (k = 0; k < K; k++)
    {
        platt_lemma_A7(err, sigma, t0, h, k, A, prec);
        _acb_vec_scalar_add_error_arb_mag(table + N*k, N, err);
    }

    arb_one(k_factorial);
    for (k = 2; k < K; k++)
    {
        row = table + N*k;
        arb_mul_ui(k_factorial, k_factorial, (ulong) k, prec);
        _acb_vec_scalar_div_arb(row, row, N, k_factorial, prec);
    }

    do_convolutions(out_table, table, S_table, N, K, prec);

    for (i = 0; i < N/2 + 1; i++)
    {
        arb_set_si(x, i);
        arb_div(x, x, B, prec);
        platt_lemma_32(err, h, t0, x, prec);
        _acb_add_error_arb_mag(out_table + i, err);
    }

    platt_lemma_B1(err, sigma, t0, h, J, prec);
    _acb_vec_scalar_add_error_arb_mag(out_table, N/2 + 1, err);

    arb_sqrt_ui(c, (ulong) J, prec);
    arb_mul_2exp_si(c, c, 1);
    arb_sub_ui(c, c, 1, prec);
    platt_lemma_B2(err, K, h, xi, prec);
    arb_mul(err, err, c, prec);
    _acb_vec_scalar_add_error_arb_mag(out_table, N/2 + 1, err);

    for (i = 1; i < N/2; i++)
    {
        acb_conj(out_table + N - i, out_table + i);
    }

    platt_lemma_A9(err, sigma, t0, h, A, prec);
    _acb_vec_scalar_add_error_arb_mag(out_table, N, err);

    acb_dft_inverse_precomp(out, out_table, pre_N, prec);
    _acb_vec_scalar_mul_arb(out, out, N, A, prec);
    for (i = 0; i < N/2; i++)
    {
        acb_swap(out + i, out + i + N/2);
    }

    platt_lemma_A11(err, t0, h, B, prec);
    _acb_vec_scalar_add_error_arb_mag(out, N, err);

    for (i = 0; i < N; i++)
    {
        n = i - N/2;
        hardy_z_per_platt_f(ratio, h, t0, n, N, B_slong, prec);
        acb_mul_arb(out + i, out + i, ratio, prec);
    }

finish:

    arb_clear(B);
    arb_clear(t);
    arb_clear(x);
    arb_clear(k_factorial);
    arb_clear(err);
    arb_clear(ratio);
    arb_clear(c);
    arb_clear(xi);
    acb_clear(z);
    _acb_vec_clear(table, K*N);
    _acb_vec_clear(S_table, K*N);
    _acb_vec_clear(out_table, N);
    acb_dft_precomp_clear(pre_N);

    return result;
}

void compare_hardy_z(slong B_slong, const arb_t A, slong N, const arb_t t0,
        const arb_t h, slong J, slong K, slong sigma, slong prec)
{
    slong i, n;
    arb_t B, t, x;
    acb_t z;
    acb_ptr observed_arr, expected_arr;
    slong digits = 12;

    arb_init(B);
    arb_init(t);
    arb_init(x);
    acb_init(z);
    observed_arr = _acb_vec_init(N);
    expected_arr = _acb_vec_init(N);

    arb_set_si(B, B_slong);

    if (!platt_multi_evaluation(observed_arr,
                B_slong, A, N, t0, h, J, K, sigma, prec))
    {
        flint_abort();
    }

    for (i = 0; i < N; i++)
    {
        n = i - N/2;
        arb_set_si(t, n);
        arb_div(t, t, A, prec);
        acb_set_arb(z, t0);
        acb_add_arb(z, z, t, prec);
        acb_dirichlet_hardy_z(expected_arr + i, z, NULL, NULL, 1, prec);
    }

    flint_printf("preliminary implementation of Platt's multi eval\n\n");
    for (i = 0; i < N; i++)
    {
        flint_printf("i: %ld\n", i);

        n = i - N/2;
        arb_set_si(t, n);
        arb_div(t, t, A, prec);
        arb_add(x, t, t0, prec);

        flint_printf("height: ");
        arb_printn(x, digits, 0);
        flint_printf("\n");

        flint_printf("multi Z : ");
        acb_printn(observed_arr + i, digits, 0);
        flint_printf("\n");

        flint_printf("single Z: ");
        acb_printn(expected_arr + i, digits, 0);
        flint_printf("\n");

        flint_printf("\n");
    }

    arb_clear(B);
    arb_clear(t);
    arb_clear(x);
    acb_clear(z);
    _acb_vec_clear(observed_arr, N);
    _acb_vec_clear(expected_arr, N);
}

int main()
{
    arb_t A, B, t0, h;
    slong K = 40;
    slong J = 10000;
    slong A_slong = 8;
    slong B_slong = 512;
    slong N = A_slong * B_slong;
    slong prec = 300;
    slong sigma = 515;

    arb_init(A);
    arb_init(B);
    arb_init(t0);
    arb_init(h);

    arb_set_ui(A, A_slong);
    arb_set_ui(B, B_slong);
    arb_set_ui(t0, 1000000);
    arb_set_ui(h, 10);

    compare_hardy_z(B_slong, A, N, t0, h, J, K, sigma, prec);

    arb_clear(A);
    arb_clear(B);
    arb_clear(t0);
    arb_clear(h);

    flint_cleanup();
    return 0;
}
