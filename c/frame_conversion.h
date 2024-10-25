// frame_conversion.h

#ifndef FRAME_CONVERSION_H
#define FRAME_CONVERSION_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct constant{
    double primary_mu; 
    double secondary_mu; 
} constant;

double dot_product(double* a, double* b, int dim);

double* cross_product(double* a, double* b);

void swap_rows(double matrix[6][6], int row1, int row2);

void swap_vector_rows(double vector[6][6], int row1, int row2);

void invert_matrix(double matrix[6][6], double inverse[6][6]);

void matrix_vector_multiply(const double matrix[6][6], const double vector[6], double result[6]); 

double* AGI_Qmethod_J2000_rv_to_RotateBarycenter_rv(double* rv_sec, double* rv_I, constant c);

double* AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(double* rv_sec, double* RV_R_barycenter_unitless, constant c);

#endif // FRAME_CONVERSION_H