/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
_acb_dirichlet_dft_precomp_init(acb_dirichlet_dft_pre_t pre, slong dv, acb_ptr z, slong dz, slong len, slong prec)
{
    if (len <= 1)
    {
        flint_printf("precomp init: trivial group. abort.\n");
        abort();
    }
    if (n_is_prime(len))
    {
        /* TODO: need convolution if len is large */
        pre->type = DFT_POL;
        _acb_dirichlet_dft_pol_init(pre->t.pol, dv, z, dz, len, prec);
    }
    else
    {
        n_factor_t fac;

        n_factor_init(&fac);
        n_factor(&fac, len, 1);

        if (fac.num == 1)
        {
            /* TODO: could be p^e, or 2^e */
            pre->type = DFT_CYC;
            _acb_dirichlet_dft_cyc_init_z_fac(pre->t.cyc, fac, dv, z, dz, len, prec);
        }
        else
        {
            pre->type = DFT_CRT;
            acb_dirichlet_dft_crt_init(pre->t.crt, len, prec);
        }
    }
}

void
acb_dirichlet_dft_precomp_init(acb_dirichlet_dft_pre_t pre, slong len, slong prec)
{
    _acb_dirichlet_dft_precomp_init(pre, 1, NULL, 0, len, prec);
}

void
acb_dirichlet_dft_precomp_clear(acb_dirichlet_dft_pre_t pre)
{
    switch (pre->type)
    {
        case DFT_POL:
            acb_dirichlet_dft_pol_clear(pre->t.pol);
            break;
        case DFT_CYC:
            acb_dirichlet_dft_cyc_clear(pre->t.cyc);
            break;
        case DFT_PROD:
            acb_dirichlet_dft_prod_clear(pre->t.prod);
            break;
        case DFT_CRT:
            acb_dirichlet_dft_crt_clear(pre->t.crt);
            break;
        default:
            flint_printf("acb_dirichlet_dft_clear: unknown strategy code %i\n", pre->type);
            abort();
            break;
    }
}

void
acb_dirichlet_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_pre_t pre, slong prec)
{
    switch (pre->type)
    {
        case DFT_POL:
            acb_dirichlet_dft_pol_precomp(w, v, pre->t.pol, prec);
            break;
        case DFT_CYC:
            acb_dirichlet_dft_cyc_precomp(w, v, pre->t.cyc, prec);
            break;
        case DFT_PROD:
            acb_dirichlet_dft_prod_precomp(w, v, pre->t.prod, prec);
            break;
        case DFT_CRT:
            acb_dirichlet_dft_crt_precomp(w, v, pre->t.crt, prec);
            break;
        default:
            flint_printf("acb_dirichlet_dft_precomp: unknown strategy code %i\n", pre->type);
            abort();
            break;
    }
}
