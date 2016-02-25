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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "fmpz_mat_extras.h"


/* fixed-capacity stack of fixed-precision signed flint integers */
typedef struct
{
    slong *data;
    slong capacity;
    slong size;
} _si_stack_struct;
typedef _si_stack_struct _si_stack_t[1];

static void
_si_stack_init(_si_stack_t S, slong capacity)
{
    S->data = flint_malloc(capacity * sizeof(slong));
    S->capacity = capacity;
    S->size = 0;
}

static void
_si_stack_clear(_si_stack_t S)
{
    flint_free(S->data);
}

static void
_si_stack_push(_si_stack_t S, slong x)
{
    if (S->size >= S->capacity) abort(); /* assert */
    S->data[S->size++] = x;
}

static slong
_si_stack_pop(_si_stack_t S)
{
    slong x;
    if (S->size <= 0) abort(); /* assert */
    x = S->data[S->size - 1];
    S->size--;
    return x;
}


/* struct for Tarjan's strongly connected components algorithm */
typedef struct
{
    slong *_index;
    slong *_lowlink;
    int *_onstack;
    _si_stack_t _S;
    slong _nsccs;
    slong _dim;
    slong _idx;
} _tarjan_struct;
typedef _tarjan_struct _tarjan_t[1];

static const slong _tarjan_UNDEFINED = -1;

static slong *
_tarjan_index(_tarjan_t t, slong i)
{
    return t->_index + i;
}

static slong *
_tarjan_lowlink(_tarjan_t t, slong i)
{
    return t->_lowlink + i;
}

static int
_tarjan_onstack(_tarjan_t t, slong i)
{
    return t->_onstack[i];
}

static void
_tarjan_push(_tarjan_t t, slong v)
{
    _si_stack_push(t->_S, v);
    t->_onstack[v] = 1;
}

static slong
_tarjan_pop(_tarjan_t t)
{
    slong v;
    v = _si_stack_pop(t->_S);
    t->_onstack[v] = 0;
    return v;
}

static slong
_tarjan_next_scc(_tarjan_t t)
{
    return t->_nsccs++;
}

static slong
_tarjan_next_idx(_tarjan_t t)
{
    return t->_idx++;
}

static void
_tarjan_init(_tarjan_t t, slong dim)
{
    slong i;
    t->_index = flint_calloc(dim, sizeof(slong));
    t->_lowlink = flint_calloc(dim, sizeof(slong));
    t->_onstack = flint_calloc(dim, sizeof(int));
    _si_stack_init(t->_S, dim);
    t->_dim = dim;
    t->_nsccs = 0;
    t->_idx = 0;
    for (i = 0; i < dim; i++)
    {
        t->_index[i] = _tarjan_UNDEFINED;
    }
}

static void
_tarjan_clear(_tarjan_t t)
{
    flint_free(t->_index);
    flint_free(t->_lowlink);
    flint_free(t->_onstack);
    _si_stack_clear(t->_S);
}

static void
_tarjan_strongconnect(slong *sccs, _tarjan_t t, const fmpz_mat_t A, slong v)
{
    slong idx, w, scc;

    idx = _tarjan_next_idx(t);
    *_tarjan_index(t, v) = idx;
    *_tarjan_lowlink(t, v) = idx;
    _tarjan_push(t, v);
    for (w = 0; w < t->_dim; w++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(A, v, w)))
        {
            if (*_tarjan_index(t, w) == _tarjan_UNDEFINED)
            {
                _tarjan_strongconnect(sccs, t, A, w);
                *_tarjan_lowlink(t, v) = FLINT_MIN(
                        *_tarjan_lowlink(t, v), *_tarjan_lowlink(t, w));
            }
            else if (_tarjan_onstack(t, w))
            {
                *_tarjan_lowlink(t, v) = FLINT_MIN(
                        *_tarjan_lowlink(t, v), *_tarjan_index(t, w));
            }
        }
    }
    if (*_tarjan_lowlink(t, v) == *_tarjan_index(t, v))
    {
        scc = _tarjan_next_scc(t);
        while (w != v)
        {
            w = _tarjan_pop(t);
            if (sccs[w] != _tarjan_UNDEFINED) abort(); /* assert */
            sccs[w] = scc;
        }
    }
}


/* Tarjan's strongly connected components algorithm */
static void
_fmpz_mat_get_sccs(slong *sccs, const fmpz_mat_t A)
{
    slong v, dim;
    _tarjan_t t;

    dim = fmpz_mat_nrows(A);

    if (dim != fmpz_mat_ncols(A))
    {
        flint_printf("_fmpz_mat_get_sccs: a square matrix is required!\n");
        abort();
    }

    _tarjan_init(t, dim);

    for (v = 0; v < dim; v++)
    {
        sccs[v] = _tarjan_UNDEFINED;
    }
    for (v = 0; v < dim; v++)
    {
        if (*_tarjan_index(t, v) == _tarjan_UNDEFINED)
        {
            _tarjan_strongconnect(sccs, t, A, v);
        }
    }

    _tarjan_clear(t);
}


/*
 * Condensation of a matrix.
 * This is the directed acyclic graph of strongly connected components.
 */
typedef struct
{
    slong n; /* number of vertices in the original graph */
    slong k; /* number of strongly connnected components (sccs) */
    fmpz_mat_t C; /* adjacency matrix of the sccs in the condensation */
    slong *partition; /* maps the vertex index to the scc index */
} _condensation_struct;

typedef _condensation_struct _condensation_t[1];

static void
_condensation_init(_condensation_t c, const fmpz_mat_t A)
{
    slong i, j, u, v;

    c->n = fmpz_mat_nrows(A);
    if (c->n != fmpz_mat_ncols(A))
    {
        flint_printf("_condensation_init: a square matrix is required!\n");
        abort();
    }

    c->partition = flint_malloc(c->n * sizeof(slong));

    _fmpz_mat_get_sccs(c->partition, A);

    /* count the strongly connected components */
    c->k = 0;
    for (i = 0; i < c->n; i++)
    {
        u = c->partition[i];
        c->k = FLINT_MAX(c->k, u + 1);
    }

    /*
     * Compute the adjacency matrix of the condensation.
     * This should be strict lower triangular, so that visiting the
     * vertices in increasing order corresponds to a postorder traversal.
     */
    fmpz_mat_init(c->C, c->k, c->k);
    fmpz_mat_zero(c->C);
    for (i = 0; i < c->n; i++)
    {
        for (j = 0; j < c->n; j++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
            {
                u = c->partition[i];
                v = c->partition[j];
                if (u != v)
                {
                    fmpz_one(fmpz_mat_entry(c->C, u, v));
                }
            }
        }
    }

    if (!fmpz_mat_is_lower_triangular(c->C) ||
        !fmpz_mat_is_hollow(c->C))
    {
        flint_printf("_condensation_init: unexpected matrix structure\n");
        fmpz_mat_print_pretty(c->C);
        abort();
    }
}

static void
_condensation_clear(_condensation_t c)
{
    fmpz_mat_clear(c->C);
    flint_free(c->partition);
}



typedef struct
{
    _condensation_t con;
    fmpz_mat_t T; /* transitive closure of condensation */
    fmpz_mat_t P; /* is there a cycle in any component on a path from u to v */
    fmpz_mat_t Q; /* longest path, if any, from u to v */
    int *scc_has_cycle;
} _connectivity_struct;
typedef _connectivity_struct _connectivity_t[1];

static void
_connectivity_clear(_connectivity_t c)
{
    fmpz_mat_clear(c->T);
    fmpz_mat_clear(c->P);
    fmpz_mat_clear(c->Q);
    flint_free(c->scc_has_cycle);
    _condensation_clear(c->con);
}

static void
_connectivity_init_scc_has_cycle(_connectivity_t c, const fmpz_mat_t A)
{
    slong n, i, u;
    slong *scc_size;

    n = fmpz_mat_nrows(A);
    c->scc_has_cycle = flint_calloc(n, sizeof(int));

    /*
     * If a vertex of the original graph has a loop,
     * then the strongly connected component to which it belongs has a cycle.
     */
    for (i = 0; i < n; i++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(A, i, i)))
        {
            u = c->con->partition[i];
            c->scc_has_cycle[u] = 1;
        }
    }

    /*
     * If a strongly connected component contains more than one vertex,
     * then that component has a cycle.
     */
    scc_size = flint_calloc(c->con->k, sizeof(slong));
    for (i = 0; i < n; i++)
    {
        u = c->con->partition[i];
        scc_size[u]++;
    }
    for (u = 0; u < c->con->k; u++)
    {
        if (scc_size[u] > 1)
        {
            c->scc_has_cycle[u] = 1;
        }
    }
    flint_free(scc_size);
}

static void
_connectivity_init(_connectivity_t c, const fmpz_mat_t A)
{
    slong u, v, w;
    slong k;
    slong curr, rest;

    /* compute the condensation */
    _condensation_init(c->con, A);
    k = c->con->k;

    /* check whether each scc contains a cycle */
    _connectivity_init_scc_has_cycle(c, A);

    /* compute the transitive closure of the condensation */
    fmpz_mat_init(c->T, k, k);
    fmpz_mat_transitive_closure(c->T, c->con->C);

    /*
     * Is there a walk from u to v that passes through a cycle-containing scc?
     * Cycles in the components u and v themselves are not considered.
     * Remember that the condensation is a directed acyclic graph.
     */
    fmpz_mat_init(c->P, k, k);
    fmpz_mat_zero(c->P);
    for (w = 0; w < k; w++)
    {
        if (c->scc_has_cycle[w])
        {
            for (u = 0; u < k; u++)
            {
                for (v = 0; v < k; v++)
                {
                    if (!fmpz_is_zero(fmpz_mat_entry(c->T, u, w)) &&
                        !fmpz_is_zero(fmpz_mat_entry(c->T, w, v)))
                    {
                        fmpz_one(fmpz_mat_entry(c->P, u, v));
                    }
                }
            }
        }
    }

    /*
     * What is the max length path from u to v in the condensation graph?
     * If u==v or if v is unreachable from u then let this be zero.
     * Remember that the condensation is a directed acyclic graph,
     * and that the components are indexed in a post-order traversal.
     */
    fmpz_mat_init(c->Q, k, k);
    fmpz_mat_zero(c->Q);
    for (u = 0; u < k; u++)
    {
        for (w = 0; w < k; w++)
        {
            if (!fmpz_is_zero(fmpz_mat_entry(c->con->C, u, w)))
            {
                curr = fmpz_get_si(fmpz_mat_entry(c->Q, u, w));
                fmpz_set_si(
                        fmpz_mat_entry(c->Q, u, w),
                        FLINT_MAX(curr, 1));
                for (v = 0; v < k; v++)
                {
                    if (!fmpz_is_zero(fmpz_mat_entry(c->T, w, v)))
                    {
                        rest = fmpz_get_si(fmpz_mat_entry(c->Q, w, v));
                        curr = fmpz_get_si(fmpz_mat_entry(c->Q, u, v));
                        fmpz_set_si(
                                fmpz_mat_entry(c->Q, u, v),
                                FLINT_MAX(curr, rest + 1));
                    }
                }
            }
        }
    }
}


static void
_connectivity_entrywise_nilpotence_degree(
        fmpz_t N, _connectivity_t c, slong i, slong j)
{
    slong u, v;
    u = c->con->partition[i];
    v = c->con->partition[j];
    if (u == v)
    {
        if (c->scc_has_cycle[u])
        {
            fmpz_set_si(N, -1);
        }
        else
        {
            fmpz_one(N);
        }
    }
    else if (fmpz_is_zero(fmpz_mat_entry(c->T, u, v)))
    {
        fmpz_zero(N);
    }
    else if (
            c->scc_has_cycle[u] ||
            c->scc_has_cycle[v] ||
            !fmpz_is_zero(fmpz_mat_entry(c->P, u, v)))
    {
        fmpz_set_si(N, -1);
    }
    else
    {
        fmpz_add_ui(N, fmpz_mat_entry(c->Q, u, v), 1);
    }
}

/*
 * Implemented only for entrywise non-negative square matrices.
 * The term 'entrywise nilpotence degree' is not standard
 * and may even be misleading, but it's supposed to mean the following:
 * For each i,j find the smallest non-negative k
 * so that (A^N)_{ij} = 0 when N >= k.
 * If there is no such k then set the entry to -1.
 * The matrix itself is nilpotent if and only if each entry
 * of the entrywise nilpotence degree matrix is non-negative.
 * If the matrix is nilpotent, then its nilpotence degree
 * is the maximum of the entrywise nilpotence degrees.
 */
void
fmpz_mat_entrywise_nilpotence_degree(fmpz_mat_t dest, const fmpz_mat_t src)
{
    slong i, j;
    _connectivity_t c;

    if (!fmpz_mat_is_square(src) || !fmpz_mat_is_nonnegative(src))
    {
        flint_printf("fmpz_mat_entrywise_nilpotence_degree: "
                     "a non-negative square matrix is required!\n");
        abort();
    }

    _connectivity_init(c, src);
    for (i = 0; i < c->con->n; i++)
    {
        for (j = 0; j < c->con->n; j++)
        {
            _connectivity_entrywise_nilpotence_degree(
                    fmpz_mat_entry(dest, i, j), c, i, j);
        }
    }
    _connectivity_clear(c);
}
