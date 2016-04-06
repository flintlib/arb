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

    Copyright (C) 2016 Jonathan Bober

******************************************************************************/

#include "arb_mat.h"

void arb_set_exact(arb_t x) {
    mag_zero(arb_radref(x));
}

void arf_twobytwo_diag(arf_t u1, arf_t u2, const arf_t a, const arf_t b, const arf_t d, slong prec) {
    /*
    // Compute the orthogonal matrix that diagonalizes
    //
    //    A = [a b]
    //        [b d]
    //
    // This matrix will have the form
    //
    //    U = [cos x , -sin x]
    //        [sin x, cos x]
    //
    // where the diagonal matrix is U^t A U.
    // We set u1 = cos x, u2 = -sin x.
    */
    arf_t x;

    if(arf_is_zero(b)) {
        arf_set_ui(u1, 1);
        arf_set_ui(u2, 0);
        return;
    }
    arf_init(x);

    arf_mul(u1, b, b, prec, ARF_RND_NEAR);            /* u1 = b^2 */
    arf_sub(u2, a, d, prec, ARF_RND_NEAR);            /* u2 = a - d */
    arf_mul_2exp_si(u2, u2, -1);                      /* u2 = (a - d)/2 */
    arf_mul(u2, u2, u2, prec, ARF_RND_NEAR);          /* u2 = ( (a - d)/2 )^2 */
    arf_add(u1, u1, u2, prec, ARF_RND_NEAR);          /* u1 = b^2 + ( (a-d)/2 )^2 */
    arf_sqrt(u1, u1, prec, ARF_RND_NEAR);             /* u1 = sqrt(above) */

    arf_mul_2exp_si(u1, u1, 1);                       /* u1 = 2 (sqrt (above) ) */
    arf_add(u1, u1, d, prec, ARF_RND_NEAR);           /* u1 += d */
    arf_sub(u1, u1, a, prec, ARF_RND_NEAR);           /* u1 -= a */
    arf_mul_2exp_si(u1, u1, -1);                      /* u1 = (d - a)/2 + sqrt(b^2 + ( (a-d)/2 )^2) */

    arf_mul(x, u1, u1, prec, ARF_RND_NEAR);
    arf_addmul(x, b, b, prec, ARF_RND_NEAR);          /* x = u1^2 + b^2 */
    arf_sqrt(x, x, prec, ARF_RND_NEAR);               /* x = sqrt(u1^2 + b^2) */
    arf_div(u2, u1, x, prec, ARF_RND_NEAR);
    arf_div(u1, b, x, prec, ARF_RND_NEAR);
    arf_neg(u1, u1);

    arf_clear(x);
}


int _arb_mat_jacobi(arb_mat_t D, arb_mat_t P, const arb_mat_t A, slong prec) {
    /*
    // Given a d x d real symmetric matrix A, compute an orthogonal matrix
    // P and a diagonal D such that A = P D P^t = P D P^(-1).
    //
    // D should have already been initialized as a d x 1 matrix, and Pp
    // should have already been initialized as a d x d matrix.
    //
    // If the eigenvalues can be certified as unique, then a nonzero int is
    // returned, and the eigenvectors should have reasonable error bounds. If
    // the eigenvalues cannot be certified as unique, then some of the
    // eigenvectors will have infinite error radius.
    */

#define B(i,j) arb_mat_entry(B, i, j)
#define D(i) arb_mat_entry(D, i, 0)
#define P(i,j) arb_mat_entry(P, i, j)

    arf_t *B1, *B2, *row_max;
    int *row_max_indices;
    int dim, i, j, k, l;
    int iter, finished;
    arb_mat_t B;
    arf_t x1, x2;
    arf_t Gii, Gij, Gji, Gjj;

    dim = arb_mat_nrows(A);
    if(dim == 1) {
        arb_mat_set(D, A);
        arb_mat_one(P);
        return 0;
    }
    arb_mat_init(B, dim, dim);

    B1 = (arf_t*)malloc(dim * sizeof(arf_t));
    B2 = (arf_t*)malloc(dim * sizeof(arf_t));
    row_max = (arf_t*)malloc((dim - 1) * sizeof(arf_t));
    row_max_indices = (int*)malloc((dim - 1) * sizeof(int));

    for(k = 0; k < dim; k++) {
        arf_init(B1[k]);
        arf_init(B2[k]);
    }
    for(k = 0; k < dim - 1; k++) {
        arf_init(row_max[k]);
    }

    arf_init(x1);
    arf_init(x2);

    arf_init(Gii);
    arf_init(Gij);
    arf_init(Gji);
    arf_init(Gjj);

    arb_mat_set(B, A);
    arb_mat_one(P);

    for(i = 0; i < dim - 1; i++) {
        for(j = i + 1; j < dim; j++) {
            arf_abs(x1, arb_midref(B(i,j)));
            if(arf_cmp(row_max[i], x1) < 0) {
                arf_set(row_max[i], x1);
                row_max_indices[i] = j;
            }
        }
    }

    finished = 0;
    iter = 0;
    while (!finished) {
        slong bound;

        iter++;
        /* flint_printf("jacobi iter %wd\n", iter); */

        arf_zero(x1);
        i = 0;
        j = 0;
        for(k = 0; k < dim - 1; k++) {
            if(arf_cmp(x1, row_max[k]) < 0) {
                arf_set(x1, row_max[k]);
                i = k;
            }
        }
        j = row_max_indices[i];

        bound = arf_abs_bound_lt_2exp_si(x1);
        if(bound < -prec * .9) {
            finished = 1;
            break;
        }

        arf_twobytwo_diag(Gii, Gij, arb_midref(B(i,i)), arb_midref(B(i,j)), arb_midref(B(j,j)), 2*prec);
        arf_neg(Gji, Gij);
        arf_set(Gjj, Gii);

        if(arf_is_zero(Gij)) {  /* If this happens, we're */
            finished = 1;       /* not going to do any better */
            break;              /* without increasing the precision. */
        }

        for(k = 0; k < dim; k++) {
            arf_mul(B1[k], Gii, arb_midref(B(i,k)), prec, ARF_RND_NEAR);
            arf_addmul(B1[k], Gji, arb_midref(B(j,k)), prec, ARF_RND_NEAR);

            arf_mul(B2[k], Gij, arb_midref(B(i,k)), prec, ARF_RND_NEAR);
            arf_addmul(B2[k], Gjj, arb_midref(B(j,k)), prec, ARF_RND_NEAR);
        }
        for(k = 0; k < dim; k++) {
            arf_set(arb_midref(B(i,k)), B1[k]);
            arf_set(arb_midref(B(j,k)), B2[k]);
        }
        for(k = 0; k < dim; k++) {
            arf_mul(B1[k], Gii, arb_midref(B(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B1[k], Gji, arb_midref(B(k,j)), prec, ARF_RND_NEAR);

            arf_mul(B2[k], Gij, arb_midref(B(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B2[k], Gjj, arb_midref(B(k,j)), prec, ARF_RND_NEAR);
        }
        for(k = 0; k < dim; k++) {
            arf_set(arb_midref(B(k,i)), B1[k]);
            arf_set(arb_midref(B(k,j)), B2[k]);
        }

        for(k = 0; k < dim; k++) {
            arf_mul(B1[k], Gii, arb_midref(P(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B1[k], Gji, arb_midref(P(k,j)), prec, ARF_RND_NEAR);

            arf_mul(B2[k], Gij, arb_midref(P(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B2[k], Gjj, arb_midref(P(k,j)), prec, ARF_RND_NEAR);
        }
        for(k = 0; k < dim; k++) {
            arf_set(arb_midref(P(k,i)), B1[k]);
            arf_set(arb_midref(P(k,j)), B2[k]);
        }

        if(i < dim - 1)
            arf_set_ui(row_max[i], 0);
        if(j < dim - 1)
            arf_set_ui(row_max[j], 0);

        /* Update the max in any row where the maximum
        // was in a column that changed. */
        for(k = 0; k < dim - 1; k++) {
            if(row_max_indices[k] == j || row_max_indices[k] == i) {
                arf_abs(row_max[k], arb_midref(B(k,k+1)));
                row_max_indices[k] = k+1;
                for(l = k+2; l < dim; l++) {
                    arf_abs(x1, arb_midref(B(k,l)));
                    if(arf_cmp(row_max[k], x1) < 0) {
                        arf_set(row_max[k], x1);
                        row_max_indices[k] = l;
                    }
                }
            }
        }

        /* Update the max in the ith row. */
        for(k = i + 1; k < dim; k++) {
            arf_abs(x1, arb_midref(B(i, k)));
            if(arf_cmp(row_max[i], x1) < 0) {
                arf_set(row_max[i], x1);
                row_max_indices[i] = k;
            }
        }

        /* Update the max in the jth row. */
        for(k = j + 1; k < dim; k++) {
            arf_abs(x1, arb_midref(B(j, k)));
            if(arf_cmp(row_max[j], x1) < 0) {
                arf_set(row_max[j], x1);
                row_max_indices[j] = k;
            }
        }

        /*
        // Go through column i to see if any of
        // the new entries are larger than the
        // max of their row.
        // */
        for(k = 0; k < i; k++) {
            if(k == dim) continue;
            arf_abs(x1, arb_midref(B(k, i)));
            if(arf_cmp(row_max[k], x1) < 0) {
                arf_set(row_max[k], x1);
                row_max_indices[k] = i;
            }
        }

        /* And then column j. */
        for(k = 0; k < j; k++) {
            if(k == dim) continue;
            arf_abs(x1, arb_midref(B(k, j)));
            if(arf_cmp(row_max[k], x1) < 0) {
                arf_set(row_max[k], x1);
                row_max_indices[k] = j;
            }
        }
    }

    for(k = 0; k < dim; k++) {
        arb_set(D(k), B(k,k));
        arb_set_exact(D(k));
    }

    /*
    // At this point we've done that diagonalization and all that remains is
    // to certify the correctness and compute error bounds.
    // */

    {
    arb_mat_t e;
    arb_struct *error_norms;
    arb_t z1, z2;
    int unique_eigenvalues = 1;

    error_norms = _arb_vec_init(dim);

    arb_mat_init(e, dim, 1);

    arb_init(z1);
    arb_init(z2);
    for(j = 0; j < dim; j++) {
        arb_mat_set(B, A);
        for(k = 0; k < dim; k++) {
            arb_sub(B(k, k), B(k, k), D(j), prec);
        }
        for(k = 0; k < dim; k++) {
            arb_set(arb_mat_entry(e, k, 0), P(k, j));
        }
        arb_mat_frobenius_norm(z2, e, prec);
        arb_mat_mul(e, B, e, prec);
        arb_mat_frobenius_norm(error_norms+j, e, prec);

        arb_div(z2, error_norms+j, z2, prec); /* and now z1 is an upper bound for the */
                                               /* error in the eigenvalue */
        arb_add_error(D(j), z2);
    }

    for(j = 0; j < dim; j++) {
        if(j == 0) {
            arb_sub(z1, D(j), D(1), prec);
        }
        else {
            arb_sub(z1, D(j), D(0), prec);
        }
        arb_get_abs_lbound_arf(x1, z1, prec);
        for(k = 1; k < dim; k++) {
            if(k == j) continue;
            arb_sub(z1, D(j), D(k), prec);
            arb_get_abs_lbound_arf(x2, z1, prec);
            if(arf_cmp(x2, x1) < 0) {
                arf_set(x1, x2);
            }
        }
        if(arf_is_zero(x1)) {
            unique_eigenvalues = 0;
        }
        arb_div_arf(z1, error_norms+j, x1, prec);
        for(k = 0; k < dim; k++) {
            arb_add_error(P(k, j), z1);
        }
    }

    arb_mat_clear(e);
    arb_clear(z1);
    arb_clear(z2);
    _arb_vec_clear(error_norms, dim);

    arf_clear(x1);
    arf_clear(x2);
    arb_mat_clear(B);
    for(k = 0; k < dim; k++) {
        arf_clear(B1[k]);
        arf_clear(B2[k]);
    }
    for(k = 0; k < dim - 1; k++) {
        arf_clear(row_max[k]);
    }
    arf_clear(Gii);
    arf_clear(Gij);
    arf_clear(Gji);
    arf_clear(Gjj);
    free(B1);
    free(B2);
    free(row_max);
    free(row_max_indices);

    if(unique_eigenvalues) return 0;
    else return 1;
#undef B
#undef D
#undef P
    }
}

int
arb_mat_symmetric_diagonalization(
        arb_mat_t D, arb_mat_t P, const arb_mat_t A, slong prec)
{
    return _arb_mat_jacobi(D, P, A, prec);
}
