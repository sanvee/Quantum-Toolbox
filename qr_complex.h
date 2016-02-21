/* linalg/qr.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_QR_COMPLEX__
#define __GSL_QR_COMPLEX__

//#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_linalg.h>

int
gsl_linalg_complex_QR_decomp (gsl_matrix_complex * A, gsl_vector_complex * tau);

int
gsl_linalg_complex_QR_lssolve (const gsl_matrix_complex * QR,
				   const gsl_vector_complex * tau,
				   const gsl_vector_complex * b,
				   gsl_vector_complex * x,
				   gsl_vector_complex * residual);

int
gsl_linalg_complex_QR_svx (const gsl_matrix_complex * QR,
			       const gsl_vector_complex * tau,
			       gsl_vector_complex * x);

int
gsl_linalg_complex_QR_QTvec (const gsl_matrix_complex * QR,
				 const gsl_vector_complex * tau,
				 gsl_vector_complex * v);

int
gsl_linalg_complex_QR_Qvec (const gsl_matrix_complex * QR,
				const gsl_vector_complex * tau,
				gsl_vector_complex * v);

int
gsl_linalg_complex_QR_unpack (const gsl_matrix_complex * QR,
				  const gsl_vector_complex * tau,
				  gsl_matrix_complex * Q,
				  gsl_matrix_complex * R);

#endif 
