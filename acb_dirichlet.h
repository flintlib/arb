/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACB_DIRICHLET_H
#define ACB_DIRICHLET_H

#ifdef ACB_DIRICHLET_INLINES_C
#define ACB_DIRICHLET_INLINE
#else
#define ACB_DIRICHLET_INLINE static __inline__
#endif

#include "acb.h"
#include "acb_poly.h"
#include "dirichlet.h"

#ifdef __cplusplus
extern "C" {
#endif

void acb_dirichlet_powsum_term(acb_ptr res, arb_t log_prev, ulong * prev,
    const acb_t s, ulong k, int integer, int critical_line, slong len, slong prec);

void acb_dirichlet_powsum_sieved(acb_ptr z, const acb_t s, ulong n, slong len, slong prec);
void acb_dirichlet_powsum_smooth(acb_ptr z, const acb_t s, ulong n, slong len, slong prec);

void acb_dirichlet_zeta_bound(mag_t res, const acb_t s);
void acb_dirichlet_zeta_rs_f_coeffs(acb_ptr c, const arb_t p, slong N, slong prec);
void acb_dirichlet_zeta_rs_d_coeffs(arb_ptr d, const arb_t sigma, slong k, slong prec);
void acb_dirichlet_zeta_rs_bound(mag_t err, const acb_t s, slong K);
void acb_dirichlet_zeta_rs_r(acb_t res, const acb_t s, slong K, slong prec);
void acb_dirichlet_zeta_rs(acb_t res, const acb_t s, slong K, slong prec);
void acb_dirichlet_zeta(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_zeta_jet_rs(acb_ptr res, const acb_t s, slong len, slong prec);
void acb_dirichlet_zeta_jet(acb_t res, const acb_t s, int deflate, slong len, slong prec);

void acb_dirichlet_hurwitz(acb_t res, const acb_t s, const acb_t a, slong prec);

void acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, slong prec);

typedef struct
{
    acb_struct s;
    mag_struct err;
    acb_ptr coeffs;
    int deflate;
    slong A;
    slong N;
    slong K;
}
acb_dirichlet_hurwitz_precomp_struct;

typedef acb_dirichlet_hurwitz_precomp_struct acb_dirichlet_hurwitz_precomp_t[1];

void acb_dirichlet_hurwitz_precomp_init(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, slong A, slong K, slong N, slong prec);
void acb_dirichlet_hurwitz_precomp_init_num(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, double num_eval, slong prec);
void acb_dirichlet_hurwitz_precomp_clear(acb_dirichlet_hurwitz_precomp_t pre);
void acb_dirichlet_hurwitz_precomp_bound(mag_t res, const acb_t s, slong A, slong K, slong N);
void acb_dirichlet_hurwitz_precomp_eval(acb_t res, const acb_dirichlet_hurwitz_precomp_t pre, ulong p, ulong q, slong prec);
void acb_dirichlet_hurwitz_precomp_choose_param(ulong * A, ulong * K, ulong * N, const acb_t s, double num_eval, slong prec);

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec);

void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_xi(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_pairing(acb_t res, const dirichlet_group_t G, ulong m, ulong n, slong prec);
void acb_dirichlet_pairing_char(acb_t res, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b, slong prec);

/* precompute roots of unity */
typedef struct
{
    ulong order;
    ulong reduced_order;
    acb_t z;
    slong size;
    slong depth;
    acb_ptr * Z;
    int use_pow;
}
acb_dirichlet_roots_struct;

typedef acb_dirichlet_roots_struct acb_dirichlet_roots_t[1];

void acb_dirichlet_roots_init(acb_dirichlet_roots_t t, ulong order, slong num, slong prec);
void acb_dirichlet_roots_clear(acb_dirichlet_roots_t t);
void acb_dirichlet_root(acb_t z, const acb_dirichlet_roots_t t, ulong n, slong prec);

void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, slong prec);
void acb_dirichlet_chi_vec(acb_ptr v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv, slong prec);

void acb_dirichlet_arb_quadratic_powers(arb_ptr v, slong nv, const arb_t x, slong prec);
void acb_dirichlet_qseries_arb(acb_t res, acb_srcptr a, const arb_t x, slong len, slong prec);
void acb_dirichlet_qseries_arb_powers_naive(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_roots_t z, slong len, slong prec);
void acb_dirichlet_qseries_arb_powers_smallorder(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_roots_t z, slong len, slong prec);

ulong acb_dirichlet_theta_length_d(ulong q, double x, slong prec);
ulong acb_dirichlet_theta_length(ulong q, const arb_t x, slong prec);
void mag_tail_kexpk2_arb(mag_t res, const arb_t a, ulong n);

void _acb_dirichlet_theta_argument_at_arb(arb_t xt, ulong q, const arb_t t, slong prec);
void acb_dirichlet_theta_arb(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t t, slong prec);
void acb_dirichlet_ui_theta_arb(acb_t res, const dirichlet_group_t G, ulong a, const arb_t t, slong prec);

void acb_dirichlet_gauss_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_order2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_root_number_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_root_number(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_si_poly_evaluate(acb_t res, slong * v, slong len, const acb_t z, slong prec);

void acb_dirichlet_jacobi_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
ulong jacobi_one_prime(ulong p, ulong e, ulong pe, ulong cond);
void acb_dirichlet_jacobi_sum_factor(acb_t res,  const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum_gauss(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, ulong b, slong prec);

void acb_dirichlet_l_euler_product(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp,
    const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);

void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp, const dirichlet_group_t G, slong prec);

void acb_dirichlet_l_jet(acb_ptr res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec);

void _acb_dirichlet_l_series(acb_ptr res, acb_srcptr s, slong slen,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    int deflate, slong len, slong prec);

void acb_dirichlet_l_series(acb_poly_t res, const acb_poly_t s,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    int deflate, slong len, slong prec);

void acb_dirichlet_hardy_theta(acb_ptr res, const acb_t t,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec);

void acb_dirichlet_hardy_z(acb_ptr res, const acb_t t,
    const dirichlet_group_t G, const dirichlet_char_t chi,
    slong len, slong prec);

void _acb_dirichlet_hardy_theta_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);
void acb_dirichlet_hardy_theta_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);
void _acb_dirichlet_hardy_z_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);
void acb_dirichlet_hardy_z_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec);

void acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec);
void acb_dirichlet_backlund_s_bound(mag_t res, const arb_t t);

/* Discrete Fourier Transform */

void acb_dirichlet_dft_index(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec);
void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec);

/* utils */

ACB_DIRICHLET_INLINE void
acb_vec_printd(acb_srcptr vec, slong len, slong digits)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_printd(vec + i, digits), flint_printf("\n");
}

#ifdef __cplusplus
}
#endif

#endif
