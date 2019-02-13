/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_calc.h"

static void
_acb_set_arf(acb_t res, const arf_t t)
{
    acb_zero(res);
    arb_set_arf(acb_realref(res), t);
}

static int
_definite_hardy_z(arb_t res, const arf_t t, slong *pprec)
{
    int msign;
    acb_t z;
    acb_init(z);
    while (1)
    {
        _acb_set_arf(z, t);
        acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, *pprec);
        msign = arb_sgn_nonzero(acb_realref(z));
        if (msign)
        {
            break;
        }
        *pprec *= 2;
    }
    acb_get_real(res, z);
    acb_clear(z);
    return msign;
}

typedef struct _zz_node_struct
{
    arf_struct t;
    arb_struct v;
    fmpz *gram;
    slong prec;
    struct _zz_node_struct *prev;
    struct _zz_node_struct *next;
}
zz_node_struct;

typedef zz_node_struct zz_node_t[1];
typedef zz_node_struct * zz_node_ptr;
typedef const zz_node_struct * zz_node_srcptr;

static int
zz_node_is_gram_node(const zz_node_t p)
{
    return p->gram != NULL;
}

static int
zz_node_sgn(const zz_node_t p)
{
    int s = arb_sgn_nonzero(&p->v);
    if (!s)
    {
        flint_printf("unexpectedly imprecise evaluation of Z(t)\n");
        flint_abort();
    }
    return s;
}

static int
zz_node_is_good_gram_node(const zz_node_t p)
{
    if (zz_node_is_gram_node(p))
    {
        int s = zz_node_sgn(p);
        if ((s > 0 && fmpz_is_even(p->gram)) ||
            (s < 0 && fmpz_is_odd(p->gram)))
        {
            return 1;
        }
    }
    return 0;
}

static void
zz_node_init(zz_node_t p)
{
    arf_init(&p->t);
    arb_init(&p->v);
    arb_indeterminate(&p->v);
    p->prec = 0;
    p->gram = NULL;
    p->prev = NULL;
    p->next = NULL;
}

static void
zz_node_clear(zz_node_t p)
{
    arf_clear(&p->t);
    arb_clear(&p->v);
    if (p->gram)
    {
        fmpz_clear(p->gram);
        flint_free(p->gram);
    }
    p->prec = 0;
    p->gram = NULL;
    p->prev = NULL;
    p->next = NULL;
}

static int
zz_node_refine(zz_node_t p)
{
    p->prec = 2*FLINT_MAX(p->prec, 8);
    return _definite_hardy_z(&p->v, &p->t, &p->prec);
}

static void
refine_one(zz_node_t a, zz_node_t b)
{
    zz_node_refine((a->prec < b->prec) ? a : b);
}

static zz_node_ptr
create_non_gram_node(const arf_t t)
{
    zz_node_ptr p = flint_malloc(sizeof(zz_node_struct));
    zz_node_init(p);
    arf_set(&p->t, t);
    p->prec = 8;
    _definite_hardy_z(&p->v, &p->t, &p->prec);
    return p;
}

static zz_node_ptr
create_gram_node(const fmpz_t n)
{
    zz_node_ptr p;
    arb_t t, v;
    acb_t z;
    slong prec = 8;

    arb_init(t);
    arb_init(v);
    acb_init(z);

    arb_indeterminate(v);
    while (arb_contains_zero(v))
    {
        prec *= 2;
        acb_dirichlet_gram_point(t, n, NULL, NULL, prec);
        acb_set_arb(z, t);
        acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, prec);
        acb_get_real(v, z);
    }

    p = flint_malloc(sizeof(zz_node_struct));
    zz_node_init(p);
    p->gram = flint_malloc(sizeof(fmpz));
    fmpz_init(p->gram);

    /* t contains g(n) and does not contain a zero of the Z function */
    fmpz_set(p->gram, n);
    arf_set(&p->t, arb_midref(t));
    arb_set(&p->v, v);
    p->prec = prec;

    arb_clear(t);
    arb_clear(v);
    acb_clear(z);

    return p;
}

static slong
count_gram_intervals(zz_node_srcptr a, zz_node_srcptr b)
{
    slong out = 0;
    if (!a || !b)
    {
        flint_printf("a and b must be non-NULL\n");
        flint_abort();
    }
    if (!zz_node_is_good_gram_node(a) || !zz_node_is_good_gram_node(b))
    {
        flint_printf("both nodes must be good Gram points\n");
        flint_abort();
    }
    else
    {
        fmpz_t m;
        fmpz_init(m);
        fmpz_sub(m, b->gram, a->gram);
        out = fmpz_get_si(m);
        fmpz_clear(m);
    }
    return out;
}

static slong
count_sign_changes(zz_node_srcptr a, zz_node_srcptr b)
{
    zz_node_srcptr p, q;
    slong n = 0;
    if (!a || !b)
    {
        flint_printf("a and b must be non-NULL\n");
        flint_abort();
    }
    p = a;
    q = a->next;
    while (p != b)
    {
        if (!q)
        {
            flint_printf("prematurely reached end of list\n");
            flint_abort();
        }
        if (zz_node_sgn(p) != zz_node_sgn(q))
        {
            n++;
        }
        p = q;
        q = q->next;
    }
    return n;
}

static zz_node_ptr
extend_to_next_good_gram_node(zz_node_t p)
{
    fmpz_t n;
    zz_node_ptr q, r;

    fmpz_init(n);

    if (!zz_node_is_gram_node(p))
    {
        flint_printf("expected to begin at a gram point\n");
        flint_abort();
    }
    if (p->next)
    {
        flint_printf("expected to extend from the end of a list\n");
        flint_abort();
    }
    fmpz_set(n, p->gram);
    q = p;
    while (1)
    {
        fmpz_add_ui(n, n, 1);
        r = create_gram_node(n);
        q->next = r;
        r->prev = q;
        q = r;
        r = NULL;
        if (zz_node_is_good_gram_node(q))
        {
            break;
        }
    }

    fmpz_clear(n);

    return q;
}

static zz_node_ptr
extend_to_prev_good_gram_node(zz_node_t p)
{
    fmpz_t n;
    zz_node_ptr q, r;

    fmpz_init(n);

    if (!zz_node_is_gram_node(p))
    {
        flint_printf("expected to begin at a gram point\n");
        flint_abort();
    }
    if (p->prev)
    {
        flint_printf("expected to extend from the start of a list\n");
        flint_abort();
    }
    fmpz_set(n, p->gram);
    q = p;
    while (1)
    {
        fmpz_sub_ui(n, n, 1);
        r = create_gram_node(n);
        q->prev = r;
        r->next = q;
        q = r;
        r = NULL;
        if (zz_node_is_good_gram_node(q))
        {
            break;
        }
    }

    fmpz_clear(n);

    return q;
}

static void
split_interval(arb_t out,
        const arf_t t1, const arb_t v1, slong sign1,
        const arf_t t2, const arb_t v2, slong sign2, slong prec)
{
    if (sign1 == sign2)
    {
        arb_t r, a, b;
        arb_init(r);
        arb_init(a);
        arb_init(b);
        arb_div(r, v2, v1, prec);
        arb_sqrt(r, r, prec);
        arb_mul_arf(a, r, t1, prec);
        arb_add_arf(a, a, t2, prec);
        arb_add_ui(b, r, 1, prec);
        arb_div(out, a, b, prec);
        arb_clear(r);
        arb_clear(a);
        arb_clear(b);
    }
    else
    {
        arb_set_arf(out, t1);
        arb_add_arf(out, out, t2, prec);
        arb_mul_2exp_si(out, out, -1);
    }
}

static void
separate_zeros(zz_node_t a, zz_node_t b, slong zn, slong limitloop)
{
    arb_t t;
    slong loopnumber = 0;
    zz_node_ptr q, r, mid_node;

    if (a == b) return;

    arb_init(t);

    while (count_sign_changes(a, b) < zn)
    {
        if (limitloop > 0 && loopnumber >= limitloop)
        {
            break;
        }
        q = a;
        r = a->next;
        while (q != b)
        {
            if (!r)
            {
                flint_printf("prematurely reached end of list\n");
                flint_abort();
            }
            while (1)
            {
                split_interval(t,
                        &q->t, &q->v, zz_node_sgn(q),
                        &r->t, &r->v, zz_node_sgn(r),
                        FLINT_MIN(q->prec, r->prec));
                if (!arb_contains_arf(t, &q->t) &&
                    !arb_contains_arf(t, &r->t))
                {
                    break;
                }
                refine_one(q, r);
            }
            mid_node = create_non_gram_node(arb_midref(t));
            q->next = mid_node;
            mid_node->prev = q;
            mid_node->next = r;
            r->prev = mid_node;
            q = r;
            r = r->next;
        }
        loopnumber++;
    }

    arb_clear(t);
}

static void
count_up(arf_interval_t r, zz_node_srcptr p, const fmpz_t n)
{
    fmpz_t N;
    fmpz_init(N);
    fmpz_add_ui(N, p->gram, 1);
    while (1)
    {
        if (!p)
        {
            flint_printf("failed to isolate zero\n");
            flint_abort();
        }
        if (zz_node_sgn(p) != zz_node_sgn(p->next))
        {
            fmpz_add_ui(N, N, 1);
            if (fmpz_equal(N, n))
            {
                arf_set(&r->a, &p->t);
                arf_set(&r->b, &p->next->t);
                break;
            }
        }
        p = p->next;
    }
    fmpz_clear(N);
}

static void
trim(zz_node_ptr *out_a, zz_node_ptr *out_b,
        zz_node_ptr a, zz_node_ptr b, slong k)
{
    slong n;
    for (n = 0; n < k; n++)
    {
        a = a->next;
        while (!zz_node_is_good_gram_node(a))
        {
            a = a->next;
        }
        b = b->prev;
        while (!zz_node_is_good_gram_node(b))
        {
            b = b->prev;
        }
    }
    *out_a = a;
    *out_b = b;
}

static void
_isolate_large_hardy_z_zero(arf_interval_t r, const fmpz_t n)
{
    fmpz_t k;
    zz_node_ptr a, b, A, B;
    slong zn, sb, cgb, variations;
    slong loopcount = 4;

    fmpz_init(k);

    fmpz_sub_ui(k, n, 2);
    a = create_gram_node(k);
    fmpz_sub_ui(k, n, 1);
    b = create_gram_node(k);
    a->next = b;
    b->prev = a;

    if (!zz_node_is_good_gram_node(b))
        b = extend_to_next_good_gram_node(b);
    if (!zz_node_is_good_gram_node(a))
        a = extend_to_prev_good_gram_node(a);

    /* Extend the search to greater heights t. */
    sb = 0;
    cgb = 0;
    while (1)
    {
        zz_node_ptr nb;
        nb = extend_to_next_good_gram_node(b);
        zn = count_gram_intervals(b, nb);
        separate_zeros(b, nb, zn, loopcount);
        if (count_sign_changes(b, nb) >= zn)
        {
            cgb++;
            if (cgb % 2 == 0 && sb < cgb / 2)
            {
                sb = cgb / 2;
                if (acb_dirichlet_turing_method_bound(nb->gram) <= sb)
                {
                    b = nb;
                    break;
                }
            }
        }
        else
        {
            cgb = 0;
        }
        b = nb;
    }

    /* Extend the search to smaller heights t. */
    cgb = 0;
    while (1)
    {
        zz_node_ptr pa;
        pa = extend_to_prev_good_gram_node(a);
        zn = count_gram_intervals(pa, a);
        separate_zeros(pa, a, zn, loopcount);
        if (count_sign_changes(pa, a) >= zn)
        {
            cgb++;
            if (cgb == sb*2)
            {
                a = pa;
                break;
            }
        }
        else
        {
            cgb = 0;
        }
        a = pa;
    }

    trim(&A, &B, a, b, sb*2);
    zn = count_gram_intervals(A, B);
    separate_zeros(A, B, zn, loopcount);
    variations = count_sign_changes(A, B);
    if (variations > zn)
    {
        flint_printf("unexpected number of sign changes\n");
        flint_abort();
    }
    else if (variations < zn)
    {
        trim(&A, &B, a, b, sb);
        zn = count_gram_intervals(A, B);
        separate_zeros(A, B, zn, 0);
        if (count_sign_changes(A, B) != zn)
        {
            flint_printf("unexpected number of sign changes\n");
            flint_abort();
        }
    }

    count_up(r, A, n);

    while (a)
    {
        b = a;
        a = a->next;
        zz_node_clear(b);
        flint_free(b);
    }

    fmpz_clear(k);
}

static void
_isolate_medium_hardy_z_zero(arf_interval_t r, const fmpz_t n)
{
    fmpz_t k;
    zz_node_ptr a, b;
    slong zn;

    fmpz_init(k);

    fmpz_sub_ui(k, n, 2);
    a = create_gram_node(k);
    fmpz_sub_ui(k, n, 1);
    b = create_gram_node(k);
    a->next = b;
    b->prev = a;

    if (!zz_node_is_good_gram_node(a))
        a = extend_to_prev_good_gram_node(a);
    if (!zz_node_is_good_gram_node(b))
        b = extend_to_next_good_gram_node(b);
    zn = count_gram_intervals(a, b);
    separate_zeros(a, b, zn, 0);
    count_up(r, a, n);

    while (a)
    {
        b = a;
        a = a->next;
        zz_node_clear(b);
        flint_free(b);
    }

    fmpz_clear(k);
}

static void
_isolate_small_hardy_z_zero(arf_interval_t r, slong n)
{
    arb_t t;
    fmpz_t k;
    slong prec = 32;

    arb_init(t);
    fmpz_init(k);

    fmpz_set_si(k, n-2);
    acb_dirichlet_gram_point(t, k, NULL, NULL, prec);
    arf_set(&r->a, arb_midref(t));

    fmpz_set_si(k, n-1);
    acb_dirichlet_gram_point(t, k, NULL, NULL, prec);
    arf_set(&r->b, arb_midref(t));

    arb_clear(t);
    fmpz_clear(k);
}

static void
_isolate_hardy_z_zero(arf_interval_t r, const fmpz_t n)
{
    if (fmpz_cmp_ui(n, 126) <= 0) /* Gram's law applies */
    {
        _isolate_small_hardy_z_zero(r, fmpz_get_si(n));
    }
    else if (fmpz_cmp_ui(n, 13999526) <= 0) /* Rosser's rule applies */
    {
        _isolate_medium_hardy_z_zero(r, n);
    }
    else
    {
        _isolate_large_hardy_z_zero(r, n);
    }
}


static int
_partition_hardy_z(arf_interval_t L, arf_interval_t R,
        const arf_interval_t block, slong *pprec)
{
    arb_t t;
    arf_t u;
    int msign;

    arb_init(t);
    arf_init(u);

    /* Compute the midpoint */
    arf_add(u, &block->a, &block->b, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(u, u, -1);

    /* Evaluate and get sign at midpoint */
    msign = 0;
    msign = _definite_hardy_z(t, u, pprec);

    /* L, R = block, split at midpoint */
    arf_set(&L->a, &block->a);
    arf_set(&R->b, &block->b);
    arf_set(&L->b, u);
    arf_set(&R->a, u);

    arb_clear(t);
    arf_clear(u);

    return msign;
}

static void
_refine_hardy_z_zero_bisect(arf_interval_t res,
        const arf_interval_t start, slong iter)
{
    int asign, bsign, msign;
    slong i, prec = 8;
    arf_interval_t t, u;
    arb_t m, v;

    arf_interval_init(t);
    arf_interval_init(u);
    arb_init(m);
    arb_init(v);

    asign = _definite_hardy_z(v, &start->a, &prec);
    bsign = _definite_hardy_z(v, &start->b, &prec);
    if (asign == bsign)
    {
        flint_printf("isolate a zero before bisecting the interval\n");
        flint_abort();
    }
    arf_interval_set(res, start);
    for (i = 0; i < iter; i++)
    {
        msign = _partition_hardy_z(t, u, res, &prec);
        if (msign == asign)
            arf_interval_swap(res, u);
        else
            arf_interval_swap(res, t);
    }

    arf_interval_clear(t);
    arf_interval_clear(u);
    arb_clear(m);
    arb_clear(v);
}

static void
_hardy_z_zero(arb_t res, const fmpz_t n, slong prec)
{
    arf_interval_t r, s;
    slong bits;

    arf_interval_init(r);
    arf_interval_init(s);

    _isolate_hardy_z_zero(r, n);

    bits = arf_bits(&r->b);
    arb_set_interval_arf(res, &r->a, &r->b, bits + 8);
    bits = arb_rel_accuracy_bits(res);

    if (bits < prec)
    {
        _refine_hardy_z_zero_bisect(s, r, prec - bits);
        arb_set_interval_arf(res, &s->a, &s->b, prec);
    }

    arb_set_round(res, res, prec);

    arf_interval_clear(r);
    arf_interval_clear(s);
}

void
acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, slong prec)
{
    fmpz_t k;
    fmpz_init(k);
    switch (fmpz_sgn(n))
    {
        case -1:
            acb_set_d(res, 0.5);
            fmpz_neg(k, n);
            _hardy_z_zero(acb_imagref(res), k, prec);
            acb_conj(res, res);
            break;
        case 1:
            acb_set_d(res, 0.5);
            _hardy_z_zero(acb_imagref(res), n, prec);
            break;
        default:
            acb_indeterminate(res);
    }
    fmpz_clear(k);
}
