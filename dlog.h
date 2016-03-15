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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#ifndef DLOG_H
#define DLOG_H

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "padic.h"

enum
{
    DLOG_MODPE, DLOG_CRT, DLOG_POWER, DLOG_BSGS, DLOG_TABLE, DLOG_23
};

typedef struct dlog_precomp_struct dlog_precomp_struct;

/* log in (1+pZ/p^eZ): compute via p-adic log */
typedef struct
{
    ulong p;
    ulong e;
    padic_ctx_t ctx;     /* padic context */
    padic_t invlog;
    padic_t x;
    fmpz_t r;
}
dlog_1modpe_struct;

typedef dlog_1modpe_struct dlog_1modpe_t[1];

/* log in (Z/p^eZ)^* */
typedef struct
{
    ulong p;
    ulong e;
    ulong pe;
    dlog_precomp_struct * modp;
    dlog_1modpe_t modpe;
}
dlog_modpe_struct;

typedef dlog_modpe_struct dlog_modpe_t[1];

/* all logs precomputed in (Z/modZ)^ast */
typedef struct
{
    ulong mod;
    ulong * table;
}
dlog_table_struct;

typedef dlog_table_struct dlog_table_t[1];

/* bsgs table, already in flint */

typedef struct apow {
    ulong k;
    ulong ak;
} apow_t;

int apow_cmp(const apow_t * x, const apow_t * y);

typedef struct {
    nmod_t mod;
    ulong m;
    ulong am;
    ulong g;
    apow_t * table;
} dlog_bsgs_struct;

typedef dlog_bsgs_struct dlog_bsgs_t[1];
/* typedef bsgs_t dlog_bsgs_t; */

/* Pollard rho */
typedef struct {
  ulong a;
  nmod_t n;
  nmod_t mod;
  int nisprime;
} dlog_rho_struct;

typedef dlog_rho_struct dlog_rho_t[1];

/* CRT decomposition (Pohlig-Hellman) */
typedef struct
{
    nmod_t mod;
    nmod_t n;
    ulong num;
    ulong * expo;
    ulong * crt_coeffs;
    dlog_precomp_struct ** pre;
    /*
       void * pre;
       dlog_precomp_t * pre;
       */
}
dlog_crt_struct;

typedef dlog_crt_struct dlog_crt_t[1];

/* dlog when generator has prime power order */
typedef struct
{
    nmod_t mod;
    ulong p;
    ulong e;
    ulong * apk;
    dlog_precomp_struct * pre;
    /*
       void * pre;
       dlog_precomp_t * pre;
       */
}
dlog_power_struct;

typedef dlog_power_struct dlog_power_t[1];

typedef ulong dlog_order23_t[1];

/* generic decomposition */
/*typedef */
struct
dlog_precomp_struct
{
    int type;
    ulong cost;
    union
    {
        dlog_table_t table;
        dlog_bsgs_t bsgs;
        dlog_crt_t crt;
        dlog_power_t power;
        dlog_modpe_t modpe;
        dlog_order23_t order23;
    } t;
};

typedef dlog_precomp_struct dlog_precomp_t[1];

ulong dlog_order23_init(dlog_order23_t t, ulong a);
ulong dlog_table_init(dlog_table_t t, ulong a, ulong mod);
ulong dlog_crt_init(dlog_crt_t t, ulong a, ulong mod, ulong n, ulong num);
ulong dlog_power_init(dlog_power_t t, ulong a, ulong mod, ulong p, ulong e, ulong num);
ulong dlog_modpe_init(dlog_modpe_t t, ulong a, ulong p, ulong e, ulong pe, ulong num);
ulong dlog_bsgs_init(dlog_bsgs_t t, ulong a, ulong mod, ulong n, ulong m);
void dlog_1modpe_init(dlog_1modpe_t t, ulong a1, ulong p, ulong e);
void dlog_rho_init(dlog_rho_t t, ulong a, ulong mod, ulong n);
/*#define dlog_bsgs_init(t, a, n, m) bsgs_table_init(t, a, n, m)*/

void dlog_order23_clear(dlog_order23_t t);
void dlog_table_clear(dlog_table_t t);
void dlog_1modpe_clear(dlog_1modpe_t t);
void dlog_crt_clear(dlog_crt_t t);
void dlog_power_clear(dlog_power_t t);
void dlog_modpe_clear(dlog_modpe_t t);
void dlog_bsgs_clear(dlog_bsgs_t t);
void dlog_rho_clear(dlog_rho_t t);
/*#define dlog_bsgs_clear(t) bsgs_table_clear(t)*/

ulong dlog_order23(const dlog_order23_t t, ulong b);
ulong dlog_table(const dlog_table_t t, ulong b);
ulong dlog_crt(const dlog_crt_t t, ulong b);
ulong dlog_power(const dlog_power_t t, ulong b);
ulong dlog_modpe(const dlog_modpe_t t, ulong b);
ulong dlog_1modpe(const dlog_1modpe_t t, ulong b);
ulong dlog_bsgs(const dlog_bsgs_t t, ulong b);
ulong dlog_rho(const dlog_rho_t t, ulong b);
/*#define dlog_bsgs(t, b) n_discrete_log_bsgs_table(t, b)*/

void dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num);
void dlog_precomp_p_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong num);
void dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num);
void dlog_precomp_clear(dlog_precomp_t pre);
ulong dlog_precomp(const dlog_precomp_t pre, ulong b);

#define NOT_FOUND UWORD_MAX
#define LOOP_MAX_FACTOR 6
typedef struct
{
    ulong m, logm;
} log_pair_t;

void dlog_vec_loop(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_eratos_ph(ulong *v, ulong nv, ulong a, ulong va, ulong M, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_eratos(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_sieve_ph(ulong *v, ulong nv, ulong a, ulong va, ulong M, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_sieve(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_crt(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);

#endif
