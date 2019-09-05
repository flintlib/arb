/*
    Copyright (C) 2010 Juan Arias de Reyna
    Copyright (C) 2019 D.H.J. Polymath
    Copyright (C) 2019 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_calc.h"

static const slong LOOPCOUNT = 4;

/*
 * For a detailed explanation of the variant of Turing's method
 * implemented in this file, see:
 *
 * J. Arias de Reyna, "Programs for Riemann's zeta function", (J. A. J. van
 * Vonderen, Ed.) Leven met getallen : liber amicorum ter gelegenheid van de
 * pensionering van Herman te Riele, CWI (2012) 102-112,
 * https://ir.cwi.nl/pub/19724
 */

/*
 * This structure describes the local context in which Platt's scaled Lambda
 * function f may be estimated by Gaussian-windowed Whittaker-Shannon
 * interpolation.
 */
typedef struct
{
    /* grid location and shape */
    fmpz T; /* midpoint */
    slong A; /* resolution */
    slong B; /* width */

    /* interpolation tuning parameters */
    slong Ns_max; /* max number of support points per side */
    arb_struct H; /* standard deviation of Gaussian window */
    slong sigma;

    arb_ptr p; /* f evaluated at N = A*B points on the grid */
    acb_dirichlet_platt_ws_precomp_struct pre; /* precomp interpolation stuff */
}
platt_ctx_struct;

typedef platt_ctx_struct platt_ctx_t[1];
typedef platt_ctx_struct * platt_ctx_ptr;
typedef const platt_ctx_struct * platt_ctx_srcptr;

static void
platt_ctx_init(platt_ctx_t ctx,
        const fmpz_t T, slong A, slong B,
        const arb_t h, slong J, slong K, slong sigma_grid,
        slong Ns_max, const arb_t H, slong sigma_interp, slong prec)
{
    fmpz_init(&ctx->T);
    arb_init(&ctx->H);
    ctx->p = _arb_vec_init(A*B);
    ctx->A = A;
    ctx->B = B;
    ctx->Ns_max = Ns_max;
    ctx->sigma = sigma_interp;
    fmpz_set(&ctx->T, T);
    arb_set(&ctx->H, H);
    acb_dirichlet_platt_ws_precomp_init(&ctx->pre, A, H, sigma_interp, prec);
    acb_dirichlet_platt_multieval(ctx->p, T, A, B, h, J, K, sigma_grid, prec);
}

static void
platt_ctx_clear(platt_ctx_t ctx)
{
    slong N = ctx->A * ctx->B;
    fmpz_clear(&ctx->T);
    arb_clear(&ctx->H);
    _arb_vec_clear(ctx->p, N);
    acb_dirichlet_platt_ws_precomp_clear(&ctx->pre);
}

static void
platt_ctx_interpolate(arb_t res, arf_t deriv,
        const platt_ctx_t ctx, const arb_t t0, slong prec)
{
    acb_dirichlet_platt_ws_interpolation_precomp(res, deriv,
        &ctx->pre, t0, ctx->p, &ctx->T, ctx->A, ctx->B, ctx->Ns_max,
        &ctx->H, ctx->sigma, prec);
}

static void
platt_ctx_interpolate_arf(arb_t res, arf_t deriv,
        const platt_ctx_t ctx, const arf_t t0, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_arf(t, t0);
    platt_ctx_interpolate(res, deriv, ctx, t, prec);
    arb_clear(t);
}


/*
 * This structure describes a node of a doubly linked list.
 * Each node represents a height t at which Platt's scaled Lambda function
 * f(t) has been evaluated by interpolation on a grid.
 */
typedef struct _zz_node_struct
{
    arf_struct t; /* height t where v = f(t) is evaluated */
    arb_struct v; /* f(t) */
    fmpz *gram; /* Gram point index or NULL if not a Gram point */
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
        flint_printf("unexpectedly imprecise evaluation of f(t)\n");
        flint_abort();
    }
    return s;
}

/* Good Gram points are Gram points where sgn(f(g(n)))*(-1)^n > 0. */
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
    p->gram = NULL;
    p->prev = NULL;
    p->next = NULL;
}

static void
delete_list_to(zz_node_ptr head, zz_node_srcptr target)
{
    zz_node_ptr u, v;
    if (head)
    {
        if (head->prev)
        {
            flint_printf("expected the first node in the list\n");
            flint_abort();
        }
    }
    u = head;
    while (u != target)
    {
        if (u == NULL)
        {
            flint_printf("failed to find target within list\n");
            flint_abort();
        }
        v = u;
        u = u->next;
        zz_node_clear(v);
        flint_free(v);
    }
    if (u != NULL)
    {
        u->prev = NULL;
    }
}

static void
delete_list(zz_node_ptr head)
{
    delete_list_to(head, NULL);
}

/*
 * Create a node representing an evaluation of the scaled Lambda function
 * at arbitrary point t. Upon creation the sign of f(t) will
 * be known with certainty. If this is not possible then NULL is returned.
 */
static zz_node_ptr
create_non_gram_node(const arf_t t, const platt_ctx_t ctx, slong prec)
{
    zz_node_ptr p = flint_malloc(sizeof(zz_node_struct));
    zz_node_init(p);
    arf_set(&p->t, t);
    platt_ctx_interpolate_arf(&p->v, NULL, ctx, t, prec);
    if (arb_contains_zero(&p->v))
    {
        zz_node_clear(p);
        p = NULL;
    }
    return p;
}

/*
 * Create a node representing an evaluation of the scaled Lambda function
 * at the nth Gram point g(n). Upon creation a floating point number t will be
 * assigned to this node, with the property that there are no zeros of the
 * scaled Lambda function between t and the actual value of g(n).
 * The sign of f(t) will also be known with certainty, otherwise NULL
 * is returned.
 */
static zz_node_ptr
create_gram_node(const fmpz_t n, const platt_ctx_t ctx, slong prec)
{
    zz_node_ptr p = NULL;
    arb_t t, v;
    acb_t z;

    arb_init(t);
    arb_init(v);
    acb_init(z);

    acb_dirichlet_gram_point(t, n, NULL, NULL, prec + fmpz_sizeinbase(n, 2));
    acb_set_arb(z, t);
    platt_ctx_interpolate(v, NULL, ctx, t, prec);
    if (!arb_contains_zero(v))
    {
        /* t contains g(n) and does not contain a zero of the f function */
        p = flint_malloc(sizeof(zz_node_struct));
        zz_node_init(p);
        p->gram = flint_malloc(sizeof(fmpz));
        fmpz_init(p->gram);
        fmpz_set(p->gram, n);
        arf_set(&p->t, arb_midref(t));
        arb_set(&p->v, v);
    }

    arb_clear(t);
    arb_clear(v);
    acb_clear(z);
    return p;
}

/*
 * Count the number of Gram intervals between the Gram point
 * represented by node a and the Gram point represented by node b.
 * Traversing the linked list is not necessary because the Gram indices
 * of nodes a and b can be accessed directly.
 */
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

/*
 * Count the observed number of sign changes of f(t) by traversing
 * a linked list of evaluated points from node a to node b.
 */
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

/*
 * Modify a linked list that ends with node p,
 * by appending nodes representing Gram points.
 * Continue until a 'good' Gram point is found.
 * Returns nonzero on success.
 * Sets *out to the last point in the list.
 */
static int
extend_to_next_good_gram_node(
        zz_node_ptr *out, zz_node_t p, const platt_ctx_t ctx, slong prec)
{
    fmpz_t n;
    zz_node_ptr q, r;
    int result = 1;

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
        r = create_gram_node(n, ctx, prec);
        if (r)
        {
            q->next = r;
            r->prev = q;
            q = r;
            r = NULL;
            if (zz_node_is_good_gram_node(q))
            {
                break;
            }
        }
        else
        {
            result = 0;
            break;
        }
    }

    fmpz_clear(n);
    *out = q;
    return result;
}

/*
 * Modify a linked list that begins with node p,
 * by prepending nodes representing Gram points.
 * Continue until a 'good' Gram point is found.
 * Returns nonzero on success.
 * Sets *out to the first point in the list.
 */
static int
extend_to_prev_good_gram_node(zz_node_ptr *out,
        zz_node_t p, const platt_ctx_t ctx, slong prec)
{
    fmpz_t n;
    zz_node_ptr q, r;
    int result = 1;

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
        r = create_gram_node(n, ctx, prec);
        if (r)
        {
            q->prev = r;
            r->next = q;
            q = r;
            r = NULL;
            if (zz_node_is_good_gram_node(q))
            {
                break;
            }
        }
        else
        {
            result = 0;
            break;
        }
    }

    fmpz_clear(n);
    *out = q;
    return result;
}

static zz_node_ptr
_scan_to_prev_good_gram_node(zz_node_ptr p)
{
    zz_node_ptr u = p->prev;
    while (u)
    {
        if (zz_node_is_good_gram_node(u))
        {
            return u;
        }
        u = u->prev;
    }
    return NULL;
}

static zz_node_ptr
scan_to_prev_good_gram_node(zz_node_ptr p, slong count)
{
    slong i;
    zz_node_ptr u = p;
    for (i = 0; i < count; i++)
    {
        if ((u = _scan_to_prev_good_gram_node(u)) == NULL)
        {
            return NULL;
        }
    }
    return u;
}

static zz_node_ptr
_scan_to_next_good_gram_node(zz_node_ptr p)
{
    zz_node_ptr u = p->next;
    while (u)
    {
        if (zz_node_is_good_gram_node(u))
        {
            return u;
        }
        u = u->next;
    }
    return NULL;
}

static zz_node_ptr
scan_to_next_good_gram_node(zz_node_ptr p, slong count)
{
    slong i;
    zz_node_ptr u = p;
    for (i = 0; i < count; i++)
    {
        if ((u = _scan_to_next_good_gram_node(u)) == NULL)
        {
            return NULL;
        }
    }
    return u;
}

/*
 * res = (x1*w1 + x2*w2) / (w1 + w2)
 * Undefined if weights are not nonnegative.
 * If w1 and w2 are zero, the resulting interval contains x1 and x2.
 */
static void
_weighted_arithmetic_mean(arb_t res, const arf_t x1, const arf_t x2,
        const arb_t w1, const arb_t w2, slong prec)
{
    if (!arb_is_nonnegative(w1) || !arb_is_nonnegative(w2))
    {
        arb_indeterminate(res);
    }
    else if (arb_is_zero(w1) && arb_is_zero(w2))
    {
        arb_set_interval_arf(res, x1, x2, prec);
    }
    else if (arb_is_zero(w1))
    {
        arb_set_arf(res, x2);
    }
    else if (arb_is_zero(w2))
    {
        arb_set_arf(res, x1);
    }
    else if (arb_is_exact(w1) && arb_is_exact(w2))
    {
        arb_t a, b;
        arb_init(a);
        arb_init(b);
        arb_mul_arf(a, w1, x1, prec);
        arb_addmul_arf(a, w2, x2, prec);
        arb_add(b, w1, w2, prec);
        arb_div(res, a, b, prec);
        arb_clear(a);
        arb_clear(b);
    }
    else
    {
        arb_t a, b, r1, r2;
        arb_init(a);
        arb_init(b);
        arb_init(r1);
        arb_init(r2);

        arb_zero(a);
        arb_zero(b);
        arb_get_lbound_arf(arb_midref(a), w1, prec);
        arb_get_ubound_arf(arb_midref(b), w2, prec);
        _weighted_arithmetic_mean(r1, x1, x2, a, b, prec);

        arb_zero(a);
        arb_zero(b);
        arb_get_ubound_arf(arb_midref(a), w1, prec);
        arb_get_lbound_arf(arb_midref(b), w2, prec);
        _weighted_arithmetic_mean(r2, x1, x2, a, b, prec);

        arb_union(res, r1, r2, prec);

        arb_clear(a);
        arb_clear(b);
        arb_clear(r1);
        arb_clear(r2);
    }
}

/*
 * Split the interval (t1, t2) into the intervals (t1, out) and (out, t2)
 * in an attempt to increase the number of observed sign changes of f(t)
 * between endpoints.
 * v1 and v2 are the scaled Lambda function values at t1 and t2 respectively.
 * sign1 and sign2 are the signs of v1 and v2 respectively.
 */
static void
split_interval(arb_t out,
        const arf_t t1, const arb_t v1, slong sign1,
        const arf_t t2, const arb_t v2, slong sign2, slong prec)
{
    if (sign1 == sign2)
    {
        /*
         * out = (sqrt(v2/v1)*t1 + t2) / (sqrt(v2/v1) + 1)
         * We have f(t1)=v1, f(t2)=v2 where v1 and v2 have the same sign,
         * and we want to guess t between t1 and t2 so that f(t)
         * has the opposite sign. Try the vertex of a parabola that would touch
         * f(t)=0 between t1 and t2 and would pass through (t1,v1) and (t2,v2).
         */
        arb_t w1, w2;
        arb_init(w1);
        arb_init(w2);
        arb_abs(w1, v2); /* w1, v2 is deliberate */
        arb_sqrt(w1, w1, prec);
        arb_abs(w2, v1); /* w2, v1 is deliberate */
        arb_sqrt(w2, w2, prec);
        _weighted_arithmetic_mean(out, t1, t2, w1, w2, prec);
        arb_clear(w1);
        arb_clear(w2);
    }
    else
    {
        /*
         * out = (t1 + t2) / 2
         * There is already one sign change in this interval.
         * To find additional sign changes we would need to evaluate
         * at least two more points in the interval,
         * so begin by just splitting the interval in half at the midpoint.
         */
        arb_set_arf(out, t1);
        arb_add_arf(out, out, t2, prec);
        arb_mul_2exp_si(out, out, -1);
    }
}

/*
 * Add a new node between each pair of existing nodes in the linked list
 * of evaluated values of t, within the sublist demarcated by nodes a and b.
 * Returns nonzero on success.
 */
static int
intercalate(const platt_ctx_t ctx, zz_node_t a, zz_node_t b, slong prec)
{
    arb_t t;
    zz_node_ptr q, r, mid_node;
    int result = 1;

    if (a == NULL || b == NULL)
    {
        flint_printf("a and b must be non-NULL\n");
        flint_abort();
    }
    if (!zz_node_is_good_gram_node(a) || !zz_node_is_good_gram_node(b))
    {
        flint_printf("a and b must represent good Gram points\n");
        flint_abort();
    }

    if (a == b) return result;

    arb_init(t);

    q = a;
    r = a->next;
    while (q != b)
    {
        if (!r)
        {
            flint_printf("prematurely reached end of list\n");
            flint_abort();
        }
        split_interval(t,
                &q->t, &q->v, zz_node_sgn(q),
                &r->t, &r->v, zz_node_sgn(r), prec);
        if (arb_contains_arf(t, &q->t) || arb_contains_arf(t, &r->t))
        {
            result = 0;
            break;
        }
        mid_node = create_non_gram_node(arb_midref(t), ctx, prec);
        if (!mid_node)
        {
            result = 0;
            break;
        }
        q->next = mid_node;
        mid_node->prev = q;
        mid_node->next = r;
        r->prev = mid_node;
        q = r;
        r = r->next;
    }

    arb_clear(t);

    return result;
}

/*
 * Given a linked sublist beginning at U and ending at V defining function
 * evaluations at points that fully separate zeros of f(t) in the vicinity
 * of the nth zero, traverse the list until the nth zero is found.
 * Continue traversing the list until len consecutive isolating intervals
 * have been found, or until the end of the sublist is reached.
 * Return the number of isolated zeros found, starting at the nth zero.
 */
static slong
count_up_separated_zeros(arf_interval_ptr res,
        zz_node_srcptr U, zz_node_srcptr V, const fmpz_t n, slong len)
{
    if (len <= 0)
    {
        return 0;
    }
    else if (fmpz_sgn(n) < 1)
    {
        flint_printf("nonpositive indices of zeros are not supported\n");
        flint_abort();
    }
    else if (U == NULL || V == NULL)
    {
        flint_printf("U and V must not be NULL\n");
        flint_abort();
    }
    if (!zz_node_is_good_gram_node(U) || !zz_node_is_good_gram_node(V))
    {
        flint_printf("U and V must be good Gram points\n");
        flint_abort();
    }
    else
    {
        slong i = 0;
        zz_node_srcptr p = U;
        fmpz_t N, k;
        fmpz_init(N);
        fmpz_init(k);
        fmpz_add_ui(N, p->gram, 1);
        fmpz_set(k, n);
        while (p != V)
        {
            if (!p->next)
            {
                flint_printf("prematurely reached end of list\n");
                flint_abort();
            }
            if (zz_node_sgn(p) != zz_node_sgn(p->next))
            {
                fmpz_add_ui(N, N, 1);
                if (fmpz_equal(N, k))
                {
                    arf_set(&res[i].a, &p->t);
                    arf_set(&res[i].b, &p->next->t);
                    fmpz_add_ui(k, k, 1);
                    i++;
                    if (i == len)
                        break;
                }
            }
            p = p->next;
        }
        fmpz_clear(k);
        return i;
    }
    return 0;
}

/*
 * Create a small linked list defining the Gram block that is expected
 * to contain the nth zero according to the Gram heuristic.
 * Returns 0 if unable to create the Gram block.
 * The output node *p is the first node in the Gram block on success
 * or NULL on failure.
 * The output node *q is the last node in the Gram block on success
 * or NULL on failure.
 * Success does not necessarily mean that the Gram block contains the nth zero.
 */
static int
create_initial_gram_block(zz_node_ptr *p, zz_node_ptr *q,
        const platt_ctx_t ctx, const fmpz_t n, slong prec)
{
    zz_node_ptr u, v;
    fmpz_t k;
    slong result = 1;

    fmpz_init(k);
    *p = NULL;
    *q = NULL;

    fmpz_sub_ui(k, n, 2);
    u = create_gram_node(k, ctx, prec);
    if (!u)
    {
        result = 0;
        goto finish;
    }
    fmpz_sub_ui(k, n, 1);
    v = create_gram_node(k, ctx, prec);
    if (!v)
    {
        result = 0;
        goto finish;
    }
    u->next = v;
    v->prev = u;

    if (!zz_node_is_good_gram_node(u))
    {
        if (!extend_to_prev_good_gram_node(&u, u, ctx, prec))
        {
            result = 0;
            goto finish;
        }
    }
    if (!zz_node_is_good_gram_node(v))
    {
        if (!extend_to_next_good_gram_node(&v, v, ctx, prec))
        {
            result = 0;
            goto finish;
        }
    }

finish:
    if (result)
    {
        *p = u;
        *q = v;
    }
    else
    {
        delete_list(u);
    }
    return result;
}

/*
 * On failure returns 0 and output variables are NULL or zero.
 * The output variable *pu is the head of the output list.
 * The output variable *pv is the tail of the output list.
 * The first *pbound Gram blocks in the output list are certified as 'good',
 * and the list probably contains a few more trailing Gram blocks (not
 * necessarily certified as 'good').
 */
static int
create_initial_double_superblock(zz_node_ptr *pu, zz_node_ptr *pv,
        slong *pbound, const platt_ctx_t ctx, const fmpz_t n, slong prec)
{
    zz_node_ptr p, q, u, v;
    slong i, k, bound, zn;
    slong good_block_count;
    slong result = 1;

    *pu = NULL;
    *pv = NULL;
    *pbound = 0;

    if (!create_initial_gram_block(&p, &q, ctx, n, prec))
    {
        result = 0;
        goto finish;
    }

    /*
     * Add blocks in the forward direction until we have at least k Gram blocks,
     * where k is the turing method bound.
     * Note that these blocks are not necessarily good Gram blocks
     * and will not necessarily belong to the initial double superblock.
     * The bound may increase as the list is extended.
     */
    for (k = 1; k < acb_dirichlet_turing_method_bound(q->gram); k++)
    {
        if (!extend_to_next_good_gram_node(&q, q, ctx, prec))
        {
            result = 0;
            goto finish;
        }
    }
    bound = k;

    /*
     * Scan the list backwards, attempting to certify blocks as 'good'
     * and tracking the current number of consecutive good Gram blocks.
     */
    good_block_count = 0;
    v = q;
    while ((u = _scan_to_prev_good_gram_node(v)) != NULL)
    {
        zn = count_gram_intervals(u, v);
        for (i = 0; i < LOOPCOUNT && count_sign_changes(u, v) < zn; i++)
        {
            if (!intercalate(ctx, u, v, prec))
            {
                result = 0;
                goto finish;
            }
        }
        if (count_sign_changes(u, v) >= zn)
        {
            good_block_count++;
        }
        else
        {
            good_block_count = 0;
        }
        v = u;
    }

    if (v != p)
    {
        flint_printf("unexpected endpoint of backwards scan\n");
        flint_abort();
    }

    /*
     * Add blocks in the backwards direction until the number
     * of consecutive good Gram blocks is twice the computed bound.
     */
    u = p;
    v = p;
    while (good_block_count < 2*bound)
    {
        if (!extend_to_prev_good_gram_node(&u, u, ctx, prec))
        {
            result = 0;
            goto finish;
        }
        zn = count_gram_intervals(u, v);
        for (i = 0; i < LOOPCOUNT && count_sign_changes(u, v) < zn; i++)
        {
            if (!intercalate(ctx, u, v, prec))
            {
                result = 0;
                goto finish;
            }
        }
        if (count_sign_changes(u, v) >= zn)
        {
            good_block_count++;
        }
        else
        {
            good_block_count = 0;
        }
    }
    p = u;

finish:

    if (result)
    {
        *pu = p;
        *pv = q;
        *pbound = bound;
    }
    else
    {
        delete_list(p);
    }
    return result;
}


static slong
_isolate_zeros(arf_interval_ptr res,
        const platt_ctx_t ctx, const fmpz_t n, slong len, slong prec)
{
    zz_node_ptr x, y; /* Anchor nodes where N(t) is known */
    zz_node_ptr u, v; /* For certifying Gram nodes as 'good' */
    zz_node_ptr p, q;
    fmpz_t nnext;
    slong i, k, bound, zn, zc, zeros_count;

    fmpz_init(nnext);
    fmpz_set(nnext, n);

    p = NULL;
    zeros_count = 0;
    if (!create_initial_double_superblock(&p, &q, &bound, ctx, n, prec))
    {
        goto finish;
    }

    /*
     * Set the anchor to the central good Gram node of the
     * initial double superblock.
     * Delete the nodes in the list before the anchor.
     */
    x = scan_to_next_good_gram_node(p, bound);
    if (x == NULL)
    {
        flint_printf("missing or incomplete initial block\n");
        flint_abort();
    }
    delete_list_to(p, x);
    p = x;

    /*
     * Set v to the forward-most node in the list,
     * and track the number of consecutive good Gram blocks at that point.
     */
    v = scan_to_next_good_gram_node(p, bound);
    if (v == NULL)
    {
        flint_printf("missing or incomplete initial block\n");
        flint_abort();
    }
    k = 2*bound;
    u = v;
    while ((v = _scan_to_next_good_gram_node(v)) != NULL)
    {
        zn = count_gram_intervals(u, v);
        if (count_sign_changes(u, v) >= zn)
        {
            k++;
        }
        else
        {
            k = 0;
        }
        u = v;
    }
    if (u != q)
    {
        flint_printf("failed to scan the initial list\n");
        flint_abort();
    }
    v = u;

    /*
     * Iterate through Gram blocks. The central good Gram point in each
     * stretch of 2*bound consecutive good Gram blocks is an 'anchor' point
     * where the number of zeros less than that point is known.
     * Therefore the number of zeros between each pair of anchor points
     * is known. As anchor points are certified, isolate the zeros falling
     * between each pair.
     */
    while (1)
    {
        u = v;
        if (!extend_to_next_good_gram_node(&v, v, ctx, prec))
        {
            goto finish;
        }
        zn = count_gram_intervals(u, v);
        for (i = 0; i < LOOPCOUNT && count_sign_changes(u, v) < zn; i++)
        {
            if (!intercalate(ctx, u, v, prec))
            {
                goto finish;
            }
        }
        if (count_sign_changes(u, v) >= zn)
        {
            k++;
        }
        else
        {
            k = 0;
        }
        bound = acb_dirichlet_turing_method_bound(v->gram);
        if (k >= 2*bound && fmpz_cmp(x->gram, v->gram) < 0)
        {
            /* There are exactly zn zeros between the anchor points x and y. */
            y = scan_to_prev_good_gram_node(v, bound);
            if (!y)
            {
                flint_printf("failed to scan backwards to anchor point\n");
                flint_abort();
            }
            zn = count_gram_intervals(x, y);
            while (count_sign_changes(x, y) < zn)
            {
                if (!intercalate(ctx, x, y, prec))
                {
                    goto finish;
                }
            }
            zc = count_up_separated_zeros(res + zeros_count,
                    x, y, nnext, len - zeros_count);
            if (zc < 0 || zc > len - zeros_count)
            {
                flint_printf("unexpected number of isolated zeros\n");
                flint_abort();
            }
            zeros_count += zc;
            if (zeros_count == len)
            {
                goto finish;
            }
            fmpz_add_ui(nnext, nnext, zc);
            x = y;
            delete_list_to(p, x);
            p = x;
        }
    }

finish:
    fmpz_clear(nnext);
    delete_list(p);
    return zeros_count;
}


slong
_acb_dirichlet_platt_isolate_local_hardy_z_zeros(
        arf_interval_ptr res, const fmpz_t n, slong len,
        const fmpz_t T, slong A, slong B,
        const arb_t h, slong J, slong K, slong sigma_grid,
        slong Ns_max, const arb_t H, slong sigma_interp, slong prec)
{
    slong zeros_count;
    platt_ctx_t ctx;
    platt_ctx_init(ctx, T, A, B, h, J, K,
            sigma_grid, Ns_max, H, sigma_interp, prec);
    zeros_count = _isolate_zeros(res, ctx, n, len, prec);
    platt_ctx_clear(ctx);
    return zeros_count;
}


static void
_refine_local_hardy_z_zero_illinois(arb_t res,
        const platt_ctx_t ctx, const arf_t ra, const arf_t rb, slong prec)
{
    arf_t a, b, fa, fb, c, fc, t;
    arb_t z;
    slong k, nmag, abs_tol, wp;
    int asign, bsign, csign;

    arf_init(a);
    arf_init(b);
    arf_init(c);
    arf_init(fa);
    arf_init(fb);
    arf_init(fc);
    arf_init(t);
    arb_init(z);

    arf_set(a, ra);
    arf_set(b, rb);

    nmag = arf_abs_bound_lt_2exp_si(b);
    abs_tol = nmag - prec - 4;

    wp = prec + nmag + 8;
    platt_ctx_interpolate_arf(z, NULL, ctx, a, wp);
    asign = arb_sgn_nonzero(z);
    arf_set(fa, arb_midref(z));
    platt_ctx_interpolate_arf(z, NULL, ctx, b, wp);
    bsign = arb_sgn_nonzero(z);
    arf_set(fb, arb_midref(z));

    if (!asign || !bsign)
    {
        flint_printf("the function evaluations at the endpoints of the initial "
                "interval must not contain zero\n");
        flint_abort();
    }
    if (asign == bsign)
    {
        flint_printf("isolate a zero before bisecting the interval\n");
        flint_abort();
    }

    for (k = 0; k < 40; k++)
    {
        /* c = a - fa * (b - a) / (fb - fa) */
        arf_sub(c, b, a, wp, ARF_RND_NEAR);
        arf_sub(t, fb, fa, wp, ARF_RND_NEAR);
        arf_div(c, c, t, wp, ARF_RND_NEAR);
        arf_mul(c, c, fa, wp, ARF_RND_NEAR);
        arf_sub(c, a, c, wp, ARF_RND_NEAR);

        /* if c is not sandwiched between a and b,
           fall back to one bisection step */
        if (!arf_is_finite(c) ||
            !((arf_cmp(a, c) < 0 && arf_cmp(c, b) < 0) ||
              (arf_cmp(b, c) < 0 && arf_cmp(c, a) < 0)))
        {
            /* flint_printf("no sandwich (k = %wd)\n", k); */
            arf_add(c, a, b, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(c, c, -1);
        }

        platt_ctx_interpolate_arf(z, NULL, ctx, c, wp);
        csign = arb_sgn_nonzero(z);

        /* If the guess is close enough to a zero that the sign
         * cannot be determined, then use the derivative to
         * make an appropriately small interval around the guess. */
        if (!csign)
        {
            arf_t deriv, aprime, bprime, faprime, fbprime, err, delta;
            slong i, aprimesign, bprimesign;

            arf_init(deriv);
            arf_init(aprime);
            arf_init(bprime);
            arf_init(faprime);
            arf_init(fbprime);
            arf_init(err);
            arf_init(delta);

            arf_set_mag(err, arb_radref(z));
            platt_ctx_interpolate_arf(NULL, deriv, ctx, c, wp);
            arf_div(delta, err, deriv, wp, ARF_RND_NEAR);
            arf_mul_si(delta, delta, 3, wp, ARF_RND_NEAR);
            arf_mul_2exp_si(delta, delta, -1);
            arf_set(aprime, c);
            arf_set(bprime, c);

            /* When the context allows the interval endpoints to
             * be evaluated to relatively high precision,
             * this should not require more than one or two iterations. */
            for (i = 0; i < 5; i++)
            {
                arf_sub(aprime, aprime, delta, wp, ARF_RND_DOWN);
                arf_add(bprime, bprime, delta, wp, ARF_RND_UP);
                if (arf_cmp(a, b) < 0)
                {
                    if (arf_cmp(aprime, a) < 0)
                        arf_set(aprime, a);
                    if (arf_cmp(b, bprime) < 0)
                        arf_set(bprime, b);
                }
                else
                {
                    if (arf_cmp(aprime, b) < 0)
                        arf_set(aprime, b);
                    if (arf_cmp(a, bprime) < 0)
                        arf_set(bprime, a);
                }
                platt_ctx_interpolate_arf(z, NULL, ctx, aprime, wp);
                arf_set(faprime, arb_midref(z));
                aprimesign = arb_sgn_nonzero(z);
                platt_ctx_interpolate_arf(z, NULL, ctx, bprime, wp);
                arf_set(fbprime, arb_midref(z));
                bprimesign = arb_sgn_nonzero(z);
                if (aprimesign && bprimesign && aprimesign != bprimesign)
                {
                    arf_set(a, aprime);
                    arf_set(b, bprime);
                    arf_set(fa, faprime);
                    arf_set(fb, fbprime);
                    break;
                }
            }

            arf_clear(deriv);
            arf_clear(aprime);
            arf_clear(bprime);
            arf_clear(faprime);
            arf_clear(fbprime);
            arf_clear(err);
            arf_clear(delta);
            break;
        }
        arf_set(fc, arb_midref(z));

        if (csign != bsign)
        {
            arf_set(a, b);
            arf_set(fa, fb);
            asign = bsign;

            arf_set(b, c);
            arf_set(fb, fc);
            bsign = csign;
        }
        else
        {
            arf_set(b, c);
            arf_set(fb, fc);
            bsign = csign;

            arf_mul_2exp_si(fa, fa, -1);
        }

        arf_sub(t, a, b, wp, ARF_RND_DOWN);
        arf_abs(t, t);

        if (arf_cmpabs_2exp_si(t, abs_tol) < 0)
            break;
    }

    /* a and b may have changed places */
    if (arf_cmp(a, b) > 0)
        arf_swap(a, b);

    arb_set_interval_arf(res, a, b, prec);

    arf_clear(a);
    arf_clear(b);
    arf_clear(c);
    arf_clear(fa);
    arf_clear(fb);
    arf_clear(fc);
    arf_clear(t);
    arb_clear(z);
}


slong
_acb_dirichlet_platt_local_hardy_z_zeros(
        arb_ptr res, const fmpz_t n, slong len,
        const fmpz_t T, slong A, slong B,
        const arb_t h, slong J, slong K, slong sigma_grid,
        slong Ns_max, const arb_t H, slong sigma_interp, slong prec)
{
    slong zeros_count, i;
    arf_interval_ptr p;
    platt_ctx_t ctx;
    platt_ctx_init(
            ctx, T, A, B, h, J, K, sigma_grid, Ns_max, H, sigma_interp, prec);
    p = _arf_interval_vec_init(len);
    zeros_count = _isolate_zeros(p, ctx, n, len, prec);
    for (i = 0; i < zeros_count; i++)
    {
        _refine_local_hardy_z_zero_illinois(res+i, ctx, &p[i].a, &p[i].b, prec);
    }
    platt_ctx_clear(ctx);
    _arf_interval_vec_clear(p, len);
    return zeros_count;
}

static void
_arb_get_lbound_fmpz(fmpz_t z, const arb_t x, slong prec)
{
    arf_t u;
    arf_init(u);
    arb_get_lbound_arf(u, x, prec);
    arf_get_fmpz(z, u, ARF_RND_DOWN);
    arf_clear(u);
}

static void
_arb_div_si_si(arb_t res, slong a, slong b, slong prec)
{
    arb_set_si(res, a);
    arb_div_si(res, res, b, prec);
}

/* Compares f to g=a*10^b.
 * Returns a negative vlaue if f < g, positive value if g < f, otherwise 0. */
static int
_fmpz_cmp_a_10exp_b(const fmpz_t f, slong a, slong b)
{
    int result;
    fmpz_t g;
    fmpz_init(g);
    fmpz_set_ui(g, 10);
    fmpz_pow_ui(g, g, b);
    fmpz_mul_si(g, g, a);
    result = fmpz_cmp(f, g);
    fmpz_clear(g);
    return result;
}

/* Returns nonzero on success. */
static int
_heuristic_A8_J(slong *J, const fmpz_t n, slong prec)
{
    if (_fmpz_cmp_a_10exp_b(n, 1, 4) < 0 ||
        _fmpz_cmp_a_10exp_b(n, 3, 22) > 0)
    {
        return 0;
    }
    else
    {
        arb_t logn;
        double x;
        arb_init(logn);
        arb_log_fmpz(logn, n, prec);
        x = arf_get_d(arb_midref(logn), ARF_RND_NEAR);
        *J = (slong) exp(0.002133 + 0.4406*x + 0.0005188*x*x);
        arb_clear(logn);
        return 1;
    }
}

/* Returns nonzero on success. */
static int
_heuristic_A8_K(slong *K, const fmpz_t n, slong prec)
{
    if (_fmpz_cmp_a_10exp_b(n, 1, 4) < 0 ||
        _fmpz_cmp_a_10exp_b(n, 3, 22) > 0)
    {
        return 0;
    }
    else
    {
        arb_t logn;
        double x;
        arb_init(logn);
        arb_log_fmpz(logn, n, prec);
        x = arf_get_d(arb_midref(logn), ARF_RND_NEAR);
        *K = (slong) (72.92 + -0.8609*x + -0.004709*x*x);
        arb_clear(logn);
        return 1;
    }
}

/* Returns nonzero on success. */
static int
_heuristic_A8_h(arb_t h, const fmpz_t n, slong prec)
{
    if (_fmpz_cmp_a_10exp_b(n, 1, 4) < 0 ||
        _fmpz_cmp_a_10exp_b(n, 3, 22) > 0)
    {
        return 0;
    }
    else
    {
        arb_t logn;
        slong hnum;
        double x;
        arb_init(logn);
        arb_log_fmpz(logn, n, prec);
        x = arf_get_d(arb_midref(logn), ARF_RND_NEAR);
        hnum = (slong) (157.8 + 26.16*x + -1.008*x*x + 0.01542*x*x*x);
        _arb_div_si_si(h, hnum, 4, prec);
        arb_clear(logn);
        return 1;
    }
}

/* Returns nonzero on success. */
static int
_heuristic_A8_H(arb_t h, const fmpz_t n, slong prec)
{
    if (_fmpz_cmp_a_10exp_b(n, 1, 4) < 0 ||
        _fmpz_cmp_a_10exp_b(n, 3, 22) > 0)
    {
        return 0;
    }
    else
    {
        arb_t logn;
        slong Hnum;
        double x;
        arb_init(logn);
        arb_log_fmpz(logn, n, prec);
        x = arf_get_d(arb_midref(logn), ARF_RND_NEAR);
        Hnum = (slong) (28.53 + 5.828*x + -0.23386*x*x + 0.0035875*x*x*x);
        _arb_div_si_si(h, Hnum, 64, prec);
        arb_clear(logn);
        return 1;
    }
}

/* Returns nonzero on success. */
static int
_heuristic_A8_sigma_interp(slong *sigma, const fmpz_t n, slong prec)
{
    if (_fmpz_cmp_a_10exp_b(n, 1, 4) < 0 ||
        _fmpz_cmp_a_10exp_b(n, 3, 22) > 0)
    {
        return 0;
    }
    else
    {
        arb_t logn;
        double x;
        arb_init(logn);
        arb_log_fmpz(logn, n, prec);
        x = arf_get_d(arb_midref(logn), ARF_RND_NEAR);
        if (_fmpz_cmp_a_10exp_b(n, 3, 14) < 0)
        {
            *sigma = 25;
        }
        else
        {
            *sigma = (slong) (-30.47 + 2.994*x + -0.04116*x*x);
        }
        *sigma += 1 - (*sigma) % 2;
        arb_clear(logn);
        return 1;
    }
}

/* Returns nonzero on success. */
static int
_heuristic_A8_sigma_grid(slong *sigma, const fmpz_t n, slong prec)
{
    if (_fmpz_cmp_a_10exp_b(n, 1, 4) < 0 ||
        _fmpz_cmp_a_10exp_b(n, 3, 22) > 0)
    {
        return 0;
    }
    else
    {
        arb_t logn;
        double x;
        arb_init(logn);
        arb_log_fmpz(logn, n, prec);
        x = arf_get_d(arb_midref(logn), ARF_RND_NEAR);
        if (_fmpz_cmp_a_10exp_b(n, 3, 6) < 0)
        {
            *sigma = (slong) (-852.5 + 388.4*x + -13.174*x*x);
        }
        else if (_fmpz_cmp_a_10exp_b(n, 3, 18) < 0)
        {
            *sigma = (slong) (1967.5703 + 4.864*x + -0.1577*x*x);
        }
        else
        {
            *sigma = (slong) (-4010.8455 + 280.2*x + -3.335*x*x);
        }
        *sigma += 1 - (*sigma) % 2;
        arb_clear(logn);
        return 1;
    }
}


static platt_ctx_ptr
_create_heuristic_context(const fmpz_t n, slong prec)
{
    platt_ctx_ptr p = NULL;
    slong A = 8;
    slong B = 4096;
    slong Ns_max = 200;
    slong J, K;
    slong sigma_grid, sigma_interp;
    fmpz_t T, k;
    slong kbits;
    arb_t g, h, H;
    acb_struct z[2];

    fmpz_init(T);
    fmpz_init(k);
    arb_init(g);
    arb_init(h);
    arb_init(H);
    acb_init(z+0);
    acb_init(z+1);

    /* Estimate the height of the nth zero using gram points --
     * it's predicted to fall between g(n-2) and g(n-1). */
    fmpz_sub_ui(k, n, 2);
    kbits = fmpz_sizeinbase(k, 2);
    acb_dirichlet_gram_point(g, k, NULL, NULL, prec + kbits);

    /* Let T be the integer at the center of the evaluation grid.
     * The reason to add B/4 instead of B/2 is that only at most the
     * middle half of the grid points are likely to have enough precision
     * to isolate zeros. */
    _arb_get_lbound_fmpz(T, g, prec);
    fmpz_add_ui(T, T, B/4);

    if (!_heuristic_A8_J(&J, n, prec)) goto finish;
    if (!_heuristic_A8_K(&K, n, prec)) goto finish;
    if (!_heuristic_A8_sigma_interp(&sigma_interp, n, prec)) goto finish;
    if (!_heuristic_A8_sigma_grid(&sigma_grid, n, prec)) goto finish;
    if (!_heuristic_A8_h(h, n, prec)) goto finish;
    if (!_heuristic_A8_H(H, n, prec)) goto finish;

    p = malloc(sizeof(platt_ctx_struct));
    platt_ctx_init(p, T, A, B, h, J, K,
            sigma_grid, Ns_max, H, sigma_interp, prec);

finish:

    fmpz_clear(T);
    fmpz_clear(k);
    arb_clear(g);
    arb_clear(h);
    arb_clear(H);
    acb_clear(z+0);
    acb_clear(z+1);

    return p;
}


/* Returns the number of zeros found. */
slong
acb_dirichlet_platt_isolate_local_hardy_z_zeros(
        arf_interval_ptr res, const fmpz_t n, slong len, slong prec)
{
    slong zeros_count = 0;
    platt_ctx_ptr ctx = _create_heuristic_context(n, prec);
    if (ctx)
    {
        zeros_count = _isolate_zeros(res, ctx, n, len, prec);
        platt_ctx_clear(ctx);
        free(ctx);
    }
    return zeros_count;
}


/* Returns the number of zeros found. */
slong
acb_dirichlet_platt_local_hardy_z_zeros(
        arb_ptr res, const fmpz_t n, slong len, slong prec)
{
    slong zeros_count = 0;
    platt_ctx_ptr ctx;

    ctx = _create_heuristic_context(n, prec);
    if (ctx)
    {
        slong i;
        arf_interval_ptr p = _arf_interval_vec_init(len);
        zeros_count = _isolate_zeros(p, ctx, n, len, prec);
        for (i = 0; i < zeros_count; i++)
        {
            _refine_local_hardy_z_zero_illinois(
                res+i, ctx, &p[i].a, &p[i].b, prec);
        }
        _arf_interval_vec_clear(p, len);
        platt_ctx_clear(ctx);
        free(ctx);
    }
    return zeros_count;
}
