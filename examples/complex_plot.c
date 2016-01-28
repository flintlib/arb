/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_modular.h"
#include "profiler.h"

/* HLS algorithm from python's colorsys module */
static double
vv(double m1, double m2, double hue)
{
    hue = hue - floor(hue);

    if (hue < 1/6.)
        return m1 + (m2-m1)*hue*6.0;
    if (hue < 0.5)
        return m2;
    if (hue < 2/3.)
        return m1 + (m2-m1)*(2/3.-hue)*6.0;
    return m1;
}

static void
hls_to_rgb(double *R, double *G, double *B, double h, double l, double s)
{
    double m1, m2;

    if (s == 0.0)
    {
        *R = *G = *B = l;
        return;
    }

    if (l <= 0.5)
        m2 = l * (1.0+s);
    else
        m2 = l+s-(l*s);

    m1 = 2.0*l - m2;

    *R = vv(m1, m2, h + 1/3.);
    *G = vv(m1, m2, h);
    *B = vv(m1, m2, h - 1/3.);
}


#define PI 3.1415926535898

void
color_function(double * R, double * G, double * B, const acb_t z)
{
    double H, L, S;
    slong prec;
    arb_t t, u;

    if (!acb_is_finite(z) || acb_rel_accuracy_bits(z) < 4)
    {
        *R = *G = *B = 0.5;
        return;
    }

    arb_init(t);
    arb_init(u);

    prec = 32;

    arf_set_round(arb_midref(t), arb_midref(acb_realref(z)), prec, ARF_RND_DOWN);
    arf_set_round(arb_midref(u), arb_midref(acb_imagref(z)), prec, ARF_RND_DOWN);

    arb_atan2(t, u, t, prec);

    H = arf_get_d(arb_midref(t), ARF_RND_DOWN);
    H = (H + PI) / (2 * PI) + 0.5;
    H = H - floor(H);

    acb_abs(t, z, prec);

    if (arf_cmpabs_2exp_si(arb_midref(t), 200) > 0)
    {
        L = 1.0;
    }
    else if (arf_cmpabs_2exp_si(arb_midref(t), -200) < 0)
    {
        L = 0.0;
    }
    else
    {
        L = arf_get_d(arb_midref(t), ARF_RND_DOWN);
        L = 1.0 - 1.0/(1.0 + pow(L, 0.2));
    }

    S = 0.8;

    hls_to_rgb(R, G, B, H, L, S);

    arb_clear(t);
    arb_clear(u);
}

typedef void (*func_ptr)(acb_t, const acb_t, slong);

void
ai(acb_t res, const acb_t z, slong prec)
{
    acb_hypgeom_airy(res, NULL, NULL, NULL, z, prec);
}

void
bi(acb_t res, const acb_t z, slong prec)
{
    acb_hypgeom_airy(NULL, NULL, res, NULL, z, prec);
}

void
besselj(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_j(res, nu, z, prec);
    acb_clear(nu);
}

void
bessely(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_y(res, nu, z, prec);
    acb_clear(nu);
}

void
besseli(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_i(res, nu, z, prec);
    acb_clear(nu);
}

void
besselk(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_k(res, nu, z, prec);
    acb_clear(nu);
}

/* this function looks better when scaled */
void
modj(acb_t res, const acb_t z, slong prec)
{
    acb_modular_j(res, z, prec);
    acb_div_ui(res, res, 1728, prec);
}

int main(int argc, char *argv[])
{
    slong x, y, xnum, ynum, prec, i;
    double R, G, B, dxa, dxb, dya, dyb;
    FILE * fp;
    arf_t xa, xb, ya, yb;
    acb_t z, w;
    func_ptr func;

    if (argc < 2)
    {
        printf("complex_plot [-range xa xb ya yb] [-size xn yn] <func>\n\n");

        printf("Plots one of the predefined functions on [xa,xb] + [ya,yb]i\n");
        printf("using domain coloring, at a resolution of xn by yn pixels.\n\n");

        printf("Defaults parameters are [-10,10] + [-10,10]i and xn = yn = 512.\n\n");

        printf("The output is written to arbplot.ppm. If you have ImageMagick,\n");
        printf("run [convert arbplot.ppm arbplot.png] to get a PNG.\n\n");

        printf("Function codes <func> are:\n");
        printf("  gamma   - Gamma function\n");
        printf("  digamma - Digamma function\n");
        printf("  lgamma  - Logarithmic gamma function\n");
        printf("  zeta    - Riemann zeta function\n");
        printf("  erf     - Error function\n");
        printf("  ai      - Airy function Ai\n");
        printf("  bi      - Airy function Bi\n");
        printf("  besselj - Bessel function J_0\n");
        printf("  bessely - Bessel function Y_0\n");
        printf("  besseli - Bessel function I_0\n");
        printf("  besselk - Bessel function K_0\n");
        printf("  modj    - Modular j-function\n");
        printf("  modeta  - Dedekind eta function\n");
        printf("  barnesg - Barnes G-function\n");
        printf("  agm     - Arithmetic geometric mean\n\n");

        return 1;
    }

    xnum = 512;
    ynum = 512;
    dxa = dya = -10;
    dxb = dyb = 10;
    func = acb_gamma;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-size"))
        {
            xnum = atol(argv[i+1]);
            ynum = atol(argv[i+2]);
            i += 2;
        }
        else if (!strcmp(argv[i], "-range"))
        {
            dxa = atof(argv[i+1]);
            dxb = atof(argv[i+2]);
            dya = atof(argv[i+3]);
            dyb = atof(argv[i+4]);
            i += 4;
        }
        else if (!strcmp(argv[i], "gamma"))
        {
            func = acb_gamma;
        }
        else if (!strcmp(argv[i], "digamma"))
        {
            func = acb_digamma;
        }
        else if (!strcmp(argv[i], "lgamma"))
        {
            func = acb_lgamma;
        }
        else if (!strcmp(argv[i], "zeta"))
        {
            func = acb_zeta;
        }
        else if (!strcmp(argv[i], "erf"))
        {
            func = acb_hypgeom_erf;
        }
        else if (!strcmp(argv[i], "ai"))
        {
            func = ai;
        }
        else if (!strcmp(argv[i], "bi"))
        {
            func = bi;
        }
        else if (!strcmp(argv[i], "besselj"))
        {
            func = besselj;
        }
        else if (!strcmp(argv[i], "bessely"))
        {
            func = bessely;
        }
        else if (!strcmp(argv[i], "besseli"))
        {
            func = besseli;
        }
        else if (!strcmp(argv[i], "besselk"))
        {
            func = besselk;
        }
        else if (!strcmp(argv[i], "modj"))
        {
            func = modj;
        }
        else if (!strcmp(argv[i], "modeta"))
        {
            func = acb_modular_eta;
        }
        else if (!strcmp(argv[i], "barnesg"))
        {
            func = acb_barnes_g;
        }
        else if (!strcmp(argv[i], "agm"))
        {
            func = acb_agm1;
        }
        else
        {
            printf("unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    acb_init(z);
    acb_init(w);

    arf_init(xa);
    arf_init(xb);
    arf_init(ya);
    arf_init(yb);

    arf_set_d(xa, dxa);
    arf_set_d(xb, dxb);
    arf_set_d(ya, dya);
    arf_set_d(yb, dyb);

    fp = fopen("arbplot.ppm", "w");
    fprintf(fp, "P6\n%ld %ld 255\n", xnum, ynum);

    TIMEIT_ONCE_START

    for (y = ynum - 1; y >= 0; y--)
    {
        if (y % (ynum / 16) == 0)
            printf("row %ld\n", y);

        for (x = 0; x < xnum; x++)
        {
            for (prec = 30; prec < 500; prec *= 2)
            {
                arf_sub(arb_midref(acb_imagref(z)), yb, ya, prec, ARF_RND_DOWN);
                arf_mul_ui(arb_midref(acb_imagref(z)),
                    arb_midref(acb_imagref(z)), y, prec, ARF_RND_DOWN);
                arf_div_ui(arb_midref(acb_imagref(z)),
                    arb_midref(acb_imagref(z)), ynum - 1, prec, ARF_RND_DOWN);
                arf_add(arb_midref(acb_imagref(z)),
                    arb_midref(acb_imagref(z)), ya, prec, ARF_RND_DOWN);

                arf_sub(arb_midref(acb_realref(z)), xb, xa, prec, ARF_RND_DOWN);
                arf_mul_ui(arb_midref(acb_realref(z)),
                    arb_midref(acb_realref(z)), x, prec, ARF_RND_DOWN);
                arf_div_ui(arb_midref(acb_realref(z)),
                    arb_midref(acb_realref(z)), xnum - 1, prec, ARF_RND_DOWN);
                arf_add(arb_midref(acb_realref(z)),
                    arb_midref(acb_realref(z)), xa, prec, ARF_RND_DOWN);

                func(w, z, prec);

                if (acb_rel_accuracy_bits(w) > 4)
                    break;
            }

            color_function(&R, &G, &B, w);

            fputc(FLINT_MIN(255, floor(R * 255)), fp);
            fputc(FLINT_MIN(255, floor(G * 255)), fp);
            fputc(FLINT_MIN(255, floor(B * 255)), fp);
        }
    }

    TIMEIT_ONCE_STOP

    fclose(fp);

    arf_clear(xa);
    arf_clear(xb);
    arf_clear(ya);
    arf_clear(yb);

    acb_clear(z);
    acb_clear(w);

    flint_cleanup();
    return 0;
}

