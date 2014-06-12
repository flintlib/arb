/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "hypgeom.h"

static __inline__ void
fmpz_poly_evaluate_si(fmpz_t y, const fmpz_poly_t poly, long x)
{
    fmpz_set_si(y, x);
    fmpz_poly_evaluate_fmpz(y, poly, y);
}

static void
bsplit_recursive_fmpz(fmpz_t P, fmpz_t Q, fmpz_t B, fmpz_t T,
    const hypgeom_t hyp, long a, long b, int cont)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            fmpz_one(P);
            fmpz_one(Q);
        }
        else
        {
            fmpz_poly_evaluate_si(P, hyp->P, a);
            fmpz_poly_evaluate_si(Q, hyp->Q, a);
        }

        fmpz_poly_evaluate_si(B, hyp->B, a);
        fmpz_poly_evaluate_si(T, hyp->A, a);
        fmpz_mul(T, T, P);
    }
    else
    {
        long m;
        fmpz_t P2, Q2, B2, T2;

        m = (a + b) / 2;

        fmpz_init(P2);
        fmpz_init(Q2);
        fmpz_init(B2);
        fmpz_init(T2);

        bsplit_recursive_fmpz(P, Q, B, T, hyp, a, m, 1);
        bsplit_recursive_fmpz(P2, Q2, B2, T2, hyp, m, b, 1);

        if (fmpz_is_one(B) && fmpz_is_one(B2))
        {
            fmpz_mul(T, T, Q2);
            fmpz_addmul(T, P, T2);
        }
        else
        {
            fmpz_mul(T, T, B2);
            fmpz_mul(T, T, Q2);
            fmpz_mul(T2, T2, B);
            fmpz_addmul(T, P, T2);
        }

        fmpz_mul(B, B, B2);
        fmpz_mul(Q, Q, Q2);
        if (cont)
            fmpz_mul(P, P, P2);

        fmpz_clear(P2);
        fmpz_clear(Q2);
        fmpz_clear(B2);
        fmpz_clear(T2);
    }
}

static void
bsplit_recursive_fmprb(fmprb_t P, fmprb_t Q, fmprb_t B, fmprb_t T,
    const hypgeom_t hyp, long a, long b, int cont, long prec)
{
    if (b - a < 4)
    {
        fmpz_t PP, QQ, BB, TT;

        fmpz_init(PP);
        fmpz_init(QQ);
        fmpz_init(BB);
        fmpz_init(TT);

        bsplit_recursive_fmpz(PP, QQ, BB, TT, hyp, a, b, cont);

        fmprb_set_fmpz(P, PP);
        fmprb_set_fmpz(Q, QQ);
        fmprb_set_fmpz(B, BB);
        fmprb_set_fmpz(T, TT);

        fmpz_clear(PP);
        fmpz_clear(QQ);
        fmpz_clear(BB);
        fmpz_clear(TT);
    }
    else
    {
        long m;
        fmprb_t P2, Q2, B2, T2;

        m = (a + b) / 2;

        fmprb_init(P2);
        fmprb_init(Q2);
        fmprb_init(B2);
        fmprb_init(T2);

        bsplit_recursive_fmprb(P, Q, B, T, hyp, a, m, 1, prec);
        bsplit_recursive_fmprb(P2, Q2, B2, T2, hyp, m, b, 1, prec);

        if (fmprb_is_one(B) && fmprb_is_one(B2))
        {
            fmprb_mul(T, T, Q2, prec);
            fmprb_addmul(T, P, T2, prec);
        }
        else
        {
            fmprb_mul(T, T, B2, prec);
            fmprb_mul(T, T, Q2, prec);
            fmprb_mul(T2, T2, B, prec);
            fmprb_addmul(T, P, T2, prec);
        }

        fmprb_mul(B, B, B2, prec);
        fmprb_mul(Q, Q, Q2, prec);
        if (cont)
            fmprb_mul(P, P, P2, prec);

        fmprb_clear(P2);
        fmprb_clear(Q2);
        fmprb_clear(B2);
        fmprb_clear(T2);
    }
}

void
fmprb_hypgeom_sum(fmprb_t P, fmprb_t Q, const hypgeom_t hyp, long n, long prec)
{
    if (n < 1)
    {
        fmprb_zero(P);
        fmprb_one(Q);
    }
    else
    {
        fmprb_t B, T;
        fmprb_init(B);
        fmprb_init(T);
        bsplit_recursive_fmprb(P, Q, B, T, hyp, 0, n, 0, prec);
        if (!fmprb_is_one(B))
            fmprb_mul(Q, Q, B, prec);
        fmprb_swap(P, T);
        fmprb_clear(B);
        fmprb_clear(T);
    }
}

void
fmprb_hypgeom_infsum(fmprb_t P, fmprb_t Q, hypgeom_t hyp, long target_prec, long prec)
{
    mag_t err, z;
    long n;

    mag_init(err);
    mag_init(z);

    mag_set_fmpz(z, hyp->P->coeffs + hyp->P->length - 1);
    mag_div_fmpz(z, z, hyp->Q->coeffs + hyp->Q->length - 1);

    if (!hyp->have_precomputed)
    {
        hypgeom_precompute(hyp);
        hyp->have_precomputed = 1;
    }

    n = hypgeom_bound(err, hyp->r, hyp->boundC, hyp->boundD,
        hyp->boundK, hyp->MK, z, target_prec);

    fmprb_hypgeom_sum(P, Q, hyp, n, prec);

    if (fmpr_sgn(fmprb_midref(Q)) < 0)
    {
        fmprb_neg(P, P);
        fmprb_neg(Q, Q);
    }

    /* We have p/q = s + err i.e. (p + q*err)/q = s */
    {
        fmpr_t u, v;
        fmpr_init(u);
        fmpr_init(v);
        mag_get_fmpr(v, err);
        fmpr_add(u, fmprb_midref(Q), fmprb_radref(Q), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(u, u, v, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmprb_add_error_fmpr(P, u);
        fmpr_clear(u);
        fmpr_clear(v);
    }

    mag_clear(z);
    mag_clear(err);
}

static void
bsplit_recursive_arb(arb_t P, arb_t Q, arb_t B, arb_t T,
    const hypgeom_t hyp, long a, long b, int cont, long prec)
{
    if (b - a < 4)
    {
        fmpz_t PP, QQ, BB, TT;

        fmpz_init(PP);
        fmpz_init(QQ);
        fmpz_init(BB);
        fmpz_init(TT);

        bsplit_recursive_fmpz(PP, QQ, BB, TT, hyp, a, b, cont);

        arb_set_fmpz(P, PP);
        arb_set_fmpz(Q, QQ);
        arb_set_fmpz(B, BB);
        arb_set_fmpz(T, TT);

        fmpz_clear(PP);
        fmpz_clear(QQ);
        fmpz_clear(BB);
        fmpz_clear(TT);
    }
    else
    {
        long m;
        arb_t P2, Q2, B2, T2;

        m = (a + b) / 2;

        arb_init(P2);
        arb_init(Q2);
        arb_init(B2);
        arb_init(T2);

        bsplit_recursive_arb(P, Q, B, T, hyp, a, m, 1, prec);
        bsplit_recursive_arb(P2, Q2, B2, T2, hyp, m, b, 1, prec);

        if (arb_is_one(B) && arb_is_one(B2))
        {
            arb_mul(T, T, Q2, prec);
            arb_addmul(T, P, T2, prec);
        }
        else
        {
            arb_mul(T, T, B2, prec);
            arb_mul(T, T, Q2, prec);
            arb_mul(T2, T2, B, prec);
            arb_addmul(T, P, T2, prec);
        }

        arb_mul(B, B, B2, prec);
        arb_mul(Q, Q, Q2, prec);
        if (cont)
            arb_mul(P, P, P2, prec);

        arb_clear(P2);
        arb_clear(Q2);
        arb_clear(B2);
        arb_clear(T2);
    }
}

void
arb_hypgeom_sum(arb_t P, arb_t Q, const hypgeom_t hyp, long n, long prec)
{
    if (n < 1)
    {
        arb_zero(P);
        arb_one(Q);
    }
    else
    {
        arb_t B, T;
        arb_init(B);
        arb_init(T);
        bsplit_recursive_arb(P, Q, B, T, hyp, 0, n, 0, prec);
        if (!arb_is_one(B))
            arb_mul(Q, Q, B, prec);
        arb_swap(P, T);
        arb_clear(B);
        arb_clear(T);
    }
}

void
arb_hypgeom_infsum(arb_t P, arb_t Q, hypgeom_t hyp, long target_prec, long prec)
{
    mag_t err, z;
    long n;

    mag_init(err);
    mag_init(z);

    mag_set_fmpz(z, hyp->P->coeffs + hyp->P->length - 1);
    mag_div_fmpz(z, z, hyp->Q->coeffs + hyp->Q->length - 1);

    if (!hyp->have_precomputed)
    {
        hypgeom_precompute(hyp);
        hyp->have_precomputed = 1;
    }

    n = hypgeom_bound(err, hyp->r, hyp->boundC, hyp->boundD,
        hyp->boundK, hyp->MK, z, target_prec);

    arb_hypgeom_sum(P, Q, hyp, n, prec);

    if (arf_sgn(arb_midref(Q)) < 0)
    {
        arb_neg(P, P);
        arb_neg(Q, Q);
    }

    /* We have p/q = s + err i.e. (p + q*err)/q = s */
    {
        mag_t u;
        mag_init(u);
        arb_get_mag(u, Q);
        mag_mul(u, u, err);
        mag_add(arb_radref(P), arb_radref(P), u);
        mag_clear(u);
    }

    mag_clear(z);
    mag_clear(err);
}

