#include "frame_conversion.h"

double dot_product(double* a, double* b, int dim){
    double result = 0; 
    for(int i = 0; i < dim; i++){
        result += (a[i] * b[i]);
    }
    return result; 
}

double* cross_product(double* a, double* b){
    double* result = malloc(3 * sizeof(double));
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result; 
}

// Function to swap two rows of a matrix
void swap_rows(double matrix[6][6], int row1, int row2) {
    for (int col = 0; col < 6; ++col) {
        double temp = matrix[row1][col];
        matrix[row1][col] = matrix[row2][col];
        matrix[row2][col] = temp;
    }
}

// Function to swap two rows of a vector
void swap_vector_rows(double vector[6][6], int row1, int row2) {
    for (int col = 0; col < 6; ++col) {
        double temp = vector[row1][col];
        vector[row1][col] = vector[row2][col];
        vector[row2][col] = temp;
    }
}

// Function to invert a 6x6 matrix
void invert_matrix(double matrix[6][6], double inverse[6][6]) {
    // Initialize the inverse matrix as the identity matrix
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            inverse[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Perform Gaussian elimination
    for (int i = 0; i < 6; ++i) {
        // Find the row with the maximum element in the current column
        double max_val = fabs(matrix[i][i]);
        int max_row = i;
        for (int k = i + 1; k < 6; ++k) {
            if (fabs(matrix[k][i]) > max_val) {
                max_val = fabs(matrix[k][i]);
                max_row = k;
            }
        }

        // If the maximum value is zero, the matrix is singular
        if (fabs(max_val) < 1e-12) {
            printf("Matrix is singular and cannot be inverted.\n");
            return;
        }

        // Swap the current row with the row of the maximum element
        if (max_row != i) {
            swap_rows(matrix, i, max_row);
            swap_vector_rows(inverse, i, max_row);
        }

        // Normalize the current row
        double diagonal_val = matrix[i][i];
        for (int col = 0; col < 6; ++col) {
            matrix[i][col] /= diagonal_val;
            inverse[i][col] /= diagonal_val;
        }

        // Eliminate the current column in all other rows
        for (int k = 0; k < 6; ++k) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int col = 0; col < 6; ++col) {
                    matrix[k][col] -= factor * matrix[i][col];
                    inverse[k][col] -= factor * inverse[i][col];
                }
            }
        }
    }
    return;
}

void matrix_vector_multiply(const double matrix[6][6], const double vector[6], double result[6]) {
    for (int i = 0; i < 6; ++i) {
        result[i] = 0.0;  // Initialize each element of the result vector
        for (int j = 0; j < 6; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}


double* AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(double* rv_sec, double* RV_R_barycenter_unitless, constant c){

    // rotating frame unitless to rotating frame with unit
    double lstar  = pow(pow(rv_sec[0], 2.0) + pow(rv_sec[1], 2.0) + pow(rv_sec[2], 2.0), 0.5); // characteristic length
    double mustar = c.primary_mu + c.secondary_mu; // G*(characteristic u)
    double tstar  = pow(pow(lstar, 3.0) / mustar, 0.5); // characteristic time
    double R_R_barycenter[3] = {RV_R_barycenter_unitless[0] * lstar, RV_R_barycenter_unitless[1] * lstar, RV_R_barycenter_unitless[2] * lstar};
    double V_R_barycenter[3] = {RV_R_barycenter_unitless[3] * (lstar/tstar), RV_R_barycenter_unitless[4] * (lstar/tstar), RV_R_barycenter_unitless[5] * (lstar/tstar)};
    double RV_R_barycenter[6] = {R_R_barycenter[0], R_R_barycenter[1], R_R_barycenter[2],
                                 V_R_barycenter[0], V_R_barycenter[1], V_R_barycenter[2]};

    // take ephemeris of secondary position and calculate circular orbit velocity and create DCM
    double r_sec[3] = {rv_sec[0], rv_sec[1], rv_sec[2]}; // km
    double v_sec_ori[3] = {rv_sec[3], rv_sec[4], rv_sec[5]}; // km

    // r3bp requires circle orbit
    double norm_r_sec = pow(pow(r_sec[0], 2.0) + pow(r_sec[1], 2.0) + pow(r_sec[2], 2.0), 0.5); 
    double omega = pow(c.primary_mu / pow(norm_r_sec, 3.0), 0.5);
    double* H = cross_product(r_sec, v_sec_ori); 
    double norm_H =pow(pow(H[0], 2.0) + pow(H[1], 2.0) + pow(H[2], 2.0), 0.5); 
    double hhat[3] = {H[0] / norm_H, H[1] / norm_H, H[2] / norm_H}; 
    double omega_vector[3] = {omega * hhat[0], omega * hhat[1], omega * hhat[2]};
    double* v_sec = cross_product(omega_vector, r_sec); 

    // instanteneous rotating axes
    double xhat[3] = {r_sec[0] / norm_r_sec, r_sec[1] / norm_r_sec, r_sec[2] / norm_r_sec};
    double* r_sec_cross_v_sec = cross_product(r_sec, v_sec); 
    double norm_r_sec_cross_v_sec = pow(pow(r_sec_cross_v_sec[0], 2.0) + pow(r_sec_cross_v_sec[1], 2.0) + pow(r_sec_cross_v_sec[2], 2.0), 0.5); 
    double zhat[3] = {r_sec_cross_v_sec[0] / norm_r_sec_cross_v_sec, r_sec_cross_v_sec[1] / norm_r_sec_cross_v_sec, r_sec_cross_v_sec[2] / norm_r_sec_cross_v_sec};
    double* yhat = cross_product(zhat, xhat); 

    // position transfrom related rotating frame to inertial frame DCM (right to left)
    double R_DCM_I[3][3] = {{xhat[0], xhat[1], xhat[2]}, 
                            {yhat[0], yhat[1], yhat[2]}, 
                            {zhat[0], zhat[1], zhat[2]}};

    // another way express dxhat dyhat dzhat from AGI page benefit no need to mess with angular velocity and its direction
    double ratio_u = c.secondary_mu / (c.primary_mu + c.secondary_mu);
    double temp_a[3] = {v_sec[0] / norm_r_sec, v_sec[1] / norm_r_sec, v_sec[2] / norm_r_sec};
    double r_sec_dot_v_sec = dot_product(r_sec, v_sec, 3); 
    double temp_b[3] = {r_sec[0] * r_sec_dot_v_sec / pow(norm_r_sec, 3.0), r_sec[1] * r_sec_dot_v_sec / pow(norm_r_sec, 3.0), r_sec[2] * r_sec_dot_v_sec / pow(norm_r_sec, 3.0)};
    double dxhat[3] = {temp_a[0] - temp_b[0], temp_a[1] - temp_b[1], temp_a[2] - temp_b[2]};
    double* temp_c = cross_product(r_sec_cross_v_sec, v_sec);
    double temp_d[3] = {temp_c[0] / (norm_r_sec * norm_r_sec_cross_v_sec), temp_c[1] / (norm_r_sec * norm_r_sec_cross_v_sec), temp_c[2] / (norm_r_sec * norm_r_sec_cross_v_sec)};
    double temp_e[3] = {r_sec_cross_v_sec[0] / pow(norm_r_sec, 3.0) * norm_r_sec_cross_v_sec, r_sec_cross_v_sec[1] / pow(norm_r_sec, 3.0) * norm_r_sec_cross_v_sec, r_sec_cross_v_sec[2] / pow(norm_r_sec, 3.0) * norm_r_sec_cross_v_sec}; 
    double* temp_f = cross_product(r_sec_cross_v_sec, r_sec);
    double temp_g = dot_product(temp_e, temp_f, 3); 
    double dyhat[3] = {temp_d[0] - temp_g, temp_d[1] - temp_g, temp_d[2] - temp_g}; 
    double dzhat[3] = {0.0, 0.0, 0.0};
    double dQ_matrix[3][3] = {{dxhat[0], dxhat[1], dxhat[2]}, 
                              {dyhat[0], dyhat[1], dyhat[2]}, 
                              {dzhat[0], dzhat[1], dzhat[2]}};
    double Q[3][3] = {{R_DCM_I[0][0], R_DCM_I[0][1], R_DCM_I[0][2]}, 
                      {R_DCM_I[1][0], R_DCM_I[1][1], R_DCM_I[1][2]}, 
                      {R_DCM_I[2][0], R_DCM_I[2][1], R_DCM_I[2][2]}}; 
    double R0[3] = {ratio_u * r_sec[0], ratio_u * r_sec[1], ratio_u * r_sec[2]};
    double V0[3] = {ratio_u * v_sec[0], ratio_u * v_sec[1], ratio_u * v_sec[2]};
    double QdQQ_matrix_a[3][6] = {{Q[0][0], Q[0][1], Q[0][2], 0.0, 0.0, 0.0}, 
                                  {Q[1][0], Q[1][1], Q[1][2], 0.0, 0.0, 0.0}, 
                                  {Q[2][0], Q[2][1], Q[2][2], 0.0, 0.0, 0.0}}; 
    double QdQQ_matrix_b[6][6] = {{QdQQ_matrix_a[0][0], QdQQ_matrix_a[0][1], QdQQ_matrix_a[0][2], QdQQ_matrix_a[0][3], QdQQ_matrix_a[0][4], QdQQ_matrix_a[0][5]},
                                  {QdQQ_matrix_a[1][0], QdQQ_matrix_a[1][1], QdQQ_matrix_a[1][2], QdQQ_matrix_a[1][3], QdQQ_matrix_a[1][4], QdQQ_matrix_a[1][5]},
                                  {QdQQ_matrix_a[2][0], QdQQ_matrix_a[2][1], QdQQ_matrix_a[2][2], QdQQ_matrix_a[2][3], QdQQ_matrix_a[2][4], QdQQ_matrix_a[2][5]},
                                  {dQ_matrix[0][0],     dQ_matrix[0][1],     dQ_matrix[0][2],     Q[0][0],             Q[0][1],             Q[0][2]},
                                  {dQ_matrix[1][0],     dQ_matrix[1][1],     dQ_matrix[1][2],     Q[1][0],             Q[1][1],             Q[1][2]},
                                  {dQ_matrix[2][0],     dQ_matrix[2][1],     dQ_matrix[2][2],     Q[2][0],             Q[2][1],             Q[2][2]}};       
    
    double QdQQ_matrix_inverse[6][6];
    invert_matrix(QdQQ_matrix_b, QdQQ_matrix_inverse);

    // rotating barycenter frame -> shift barycenter to sun-centered inertial rv -> J2000 sun-centered inertial rv
    double QdQQ_matrix_inverse_dot_RV_R_barycenter[6];
    matrix_vector_multiply(QdQQ_matrix_inverse, RV_R_barycenter, QdQQ_matrix_inverse_dot_RV_R_barycenter);
    double* rv_I = malloc(6 * sizeof(double));

    rv_I[0] = QdQQ_matrix_inverse_dot_RV_R_barycenter[0] + R0[0];
    rv_I[1] = QdQQ_matrix_inverse_dot_RV_R_barycenter[1] + R0[1];
    rv_I[2] = QdQQ_matrix_inverse_dot_RV_R_barycenter[2] + R0[2];
    rv_I[3] = QdQQ_matrix_inverse_dot_RV_R_barycenter[3] + V0[0];
    rv_I[4] = QdQQ_matrix_inverse_dot_RV_R_barycenter[4] + V0[1];
    rv_I[5] = QdQQ_matrix_inverse_dot_RV_R_barycenter[5] + V0[2];

    free(H); 
    free(v_sec); 
    free(r_sec_cross_v_sec);
    free(yhat); 
    free(temp_c); 
    free(temp_f); 
    return rv_I; 
}   

double* AGI_Qmethod_J2000_rv_to_RotateBarycenter_rv(double* rv_sec, double* rv_I, constant c){

    double r_sec[3] = {rv_sec[0], rv_sec[1], rv_sec[2]};
    double v_sec_ori[3] = {rv_sec[3], rv_sec[4], rv_sec[5]};

    // r3bp requires circle orbit
    double norm_r_sec = pow(pow(r_sec[0], 2.0) + pow(r_sec[1], 2.0) + pow(r_sec[2], 2.0), 0.5); 
    double omega = pow(c.primary_mu / pow(norm_r_sec, 3.0), 0.5);
    double* H = cross_product(r_sec, v_sec_ori); 
    double norm_H =pow(pow(H[0], 2.0) + pow(H[1], 2.0) + pow(H[2], 2.0), 0.5); 
    double hhat[3] = {H[0] / norm_H, H[1] / norm_H, H[2] / norm_H}; 
    double omega_vector[3] = {omega * hhat[0], omega * hhat[1], omega * hhat[2]};
    double* v_sec = cross_product(omega_vector, r_sec); 

    // instanteneous rotating axes
    double xhat[3] = {r_sec[0] / norm_r_sec, r_sec[1] / norm_r_sec, r_sec[2] / norm_r_sec}; 
    double* r_sec_cross_v_sec = cross_product(r_sec, v_sec); 
    double norm_r_sec_cross_v_sec = pow(pow(r_sec_cross_v_sec[0], 2.0) + pow(r_sec_cross_v_sec[1], 2.0) + pow(r_sec_cross_v_sec[2], 2.0), 0.5); 
    double zhat[3] = {r_sec_cross_v_sec[0] / norm_r_sec_cross_v_sec, r_sec_cross_v_sec[1] / norm_r_sec_cross_v_sec, r_sec_cross_v_sec[2] / norm_r_sec_cross_v_sec};
    double* yhat = cross_product(zhat, xhat); 

    // position transfrom related rotating frame to inertial frame DCM (right to left)
    double R_DCM_I[3][3] = {{xhat[0], xhat[1], xhat[2]}, 
                            {yhat[0], yhat[1], yhat[2]}, 
                            {zhat[0], zhat[1], zhat[2]}};

    // another way express dxhat dyhat dzhat from AGI page benefit no need to mess with angular velocity and its direction
    double ratio_u = c.secondary_mu / (c.primary_mu + c.secondary_mu);
    double temp_a[3] = {v_sec[0] / norm_r_sec, v_sec[1] / norm_r_sec, v_sec[2] / norm_r_sec};
    double r_sec_dot_v_sec = dot_product(r_sec, v_sec, 3); 
    double temp_b[3] = {r_sec[0] * r_sec_dot_v_sec / pow(norm_r_sec, 3.0), r_sec[1] * r_sec_dot_v_sec / pow(norm_r_sec, 3.0), r_sec[2] * r_sec_dot_v_sec / pow(norm_r_sec, 3.0)};
    double dxhat[3] = {temp_a[0] - temp_b[0], temp_a[1] - temp_b[1], temp_a[2] - temp_b[2]};
    double* temp_c = cross_product(r_sec_cross_v_sec, v_sec);
    double temp_d[3] = {temp_c[0] / (norm_r_sec * norm_r_sec_cross_v_sec), temp_c[1] / (norm_r_sec * norm_r_sec_cross_v_sec), temp_c[2] / (norm_r_sec * norm_r_sec_cross_v_sec)};
    double temp_e[3] = {r_sec_cross_v_sec[0] / pow(norm_r_sec, 3.0) * norm_r_sec_cross_v_sec, r_sec_cross_v_sec[1] / pow(norm_r_sec, 3.0) * norm_r_sec_cross_v_sec, r_sec_cross_v_sec[2] / pow(norm_r_sec, 3.0) * norm_r_sec_cross_v_sec}; 
    double* temp_f = cross_product(r_sec_cross_v_sec, r_sec);
    double temp_g = dot_product(temp_e, temp_f, 3); 
    double dyhat[3] = {temp_d[0] - temp_g, temp_d[1] - temp_g, temp_d[2] - temp_g}; 
    double dzhat[3] = {0.0, 0.0, 0.0};
    double dQ_matrix[3][3] = {{dxhat[0], dxhat[1], dxhat[2]}, 
                              {dyhat[0], dyhat[1], dyhat[2]}, 
                              {dzhat[0], dzhat[1], dzhat[2]}};
    double Q[3][3] = {{R_DCM_I[0][0], R_DCM_I[0][1], R_DCM_I[0][2]}, 
                      {R_DCM_I[1][0], R_DCM_I[1][1], R_DCM_I[1][2]}, 
                      {R_DCM_I[2][0], R_DCM_I[2][1], R_DCM_I[2][2]}}; 
    double R0[3] = {ratio_u * r_sec[0], ratio_u * r_sec[1], ratio_u * r_sec[2]};
    double V0[3] = {ratio_u * v_sec[0], ratio_u * v_sec[1], ratio_u * v_sec[2]};
    double QdQQ_matrix_a[3][6] = {{Q[0][0], Q[0][1], Q[0][2], 0.0, 0.0, 0.0}, 
                                  {Q[1][0], Q[1][1], Q[1][2], 0.0, 0.0, 0.0}, 
                                  {Q[2][0], Q[2][1], Q[2][2], 0.0, 0.0, 0.0}}; 
    double QdQQ_matrix_b[6][6] = {{QdQQ_matrix_a[0][0], QdQQ_matrix_a[0][1], QdQQ_matrix_a[0][2], QdQQ_matrix_a[0][3], QdQQ_matrix_a[0][4], QdQQ_matrix_a[0][5]},
                                  {QdQQ_matrix_a[1][0], QdQQ_matrix_a[1][1], QdQQ_matrix_a[1][2], QdQQ_matrix_a[1][3], QdQQ_matrix_a[1][4], QdQQ_matrix_a[1][5]},
                                  {QdQQ_matrix_a[2][0], QdQQ_matrix_a[2][1], QdQQ_matrix_a[2][2], QdQQ_matrix_a[2][3], QdQQ_matrix_a[2][4], QdQQ_matrix_a[2][5]},
                                  {dQ_matrix[0][0],     dQ_matrix[0][1],     dQ_matrix[0][2],     Q[0][0],             Q[0][1],             Q[0][2]},
                                  {dQ_matrix[1][0],     dQ_matrix[1][1],     dQ_matrix[1][2],     Q[1][0],             Q[1][1],             Q[1][2]},
                                  {dQ_matrix[2][0],     dQ_matrix[2][1],     dQ_matrix[2][2],     Q[2][0],             Q[2][1],             Q[2][2]}};       
    
    double r_I[3] = {rv_I[0], rv_I[1], rv_I[2]};
    double v_I[3] = {rv_I[3], rv_I[4], rv_I[5]};
    double temp_h[6] = {r_I[0] - R0[0], r_I[1] - R0[1], r_I[2] - R0[2], v_I[0] - V0[0], v_I[1] - V0[1], v_I[2] - V0[2]}; 

    double RV_R_barycenter[6];
    matrix_vector_multiply(QdQQ_matrix_b, temp_h, RV_R_barycenter);
    double R_R_barycenter[3] = {RV_R_barycenter[0], RV_R_barycenter[1], RV_R_barycenter[2]};
    double V_R_barycenter[3] = {RV_R_barycenter[3], RV_R_barycenter[4], RV_R_barycenter[5]};

    double lstar  = pow(pow(rv_sec[0], 2.0) + pow(rv_sec[1], 2.0) + pow(rv_sec[2], 2.0), 0.5); // characteristic length
    double mustar = c.primary_mu + c.secondary_mu; // G*(characteristic u)
    double tstar  = pow(pow(lstar, 3.0) / mustar, 0.5); // characteristic time

    double R_R_barycenter_unitless[3] = {R_R_barycenter[0] / lstar, R_R_barycenter[1] / lstar, R_R_barycenter[2] / lstar};
    double V_R_barycenter_unitless[3] = {V_R_barycenter[0] / (lstar / tstar), V_R_barycenter[1] / (lstar / tstar), V_R_barycenter[2] / (lstar / tstar)};
    double* RV_R_barycenter_unitless = malloc(6 * sizeof(double));
    RV_R_barycenter_unitless[0] = R_R_barycenter_unitless[0]; 
    RV_R_barycenter_unitless[1] = R_R_barycenter_unitless[1]; 
    RV_R_barycenter_unitless[2] = R_R_barycenter_unitless[2]; 
    RV_R_barycenter_unitless[3] = V_R_barycenter_unitless[0]; 
    RV_R_barycenter_unitless[4] = V_R_barycenter_unitless[1]; 
    RV_R_barycenter_unitless[5] = V_R_barycenter_unitless[2]; 

    free(H); 
    free(v_sec); 
    free(r_sec_cross_v_sec);
    free(yhat); 
    free(temp_c); 
    free(temp_f); 
    return RV_R_barycenter_unitless; 
}   

int main(){
    constant c = {398600.4, 4902.8001}; 
    double rv_I[6] = {-3.643006231934826E+05, 1.795694265325507E+05, 1.374829002429234E+04, -4.224086396977384E-01, -1.135062395919662E+00, 5.566668169147170E-01};
    // printf("rv_I: [%.12f, %.12f, %.12f, %.12f, %.12f, %.12f]\n", rv_I[0], rv_I[1], rv_I[2], rv_I[3], rv_I[4], rv_I[5]);
    double rv_sec[6] = {-3.67980875E+05, 1.66471214E+05, 2.51732170E+04, -4.09612111E-01, -8.75695737E-01, -5.92758485E-02};
    // printf("rv_sec: [%.12f, %.12f, %.12f, %.12f, %.12f, %.12f]\n", rv_sec[0], rv_sec[1], rv_sec[2], rv_sec[3], rv_sec[4], rv_sec[5]);
    double* RV_R_barycenter_unitless = AGI_Qmethod_J2000_rv_to_RotateBarycenter_rv(rv_sec, rv_I, c); 
    printf("RV_R_barycenter_unitless: [%.12f, %.12f, %.12f, %.12f, %.12f, %.12f]\n", RV_R_barycenter_unitless[0], RV_R_barycenter_unitless[1], RV_R_barycenter_unitless[2], RV_R_barycenter_unitless[3], RV_R_barycenter_unitless[4], RV_R_barycenter_unitless[5]);
    double* rv_I_new = AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(rv_sec, RV_R_barycenter_unitless, c);
    printf("rv_I: [%.12f, %.12f, %.12f, %.12f, %.12f, %.12f]\n", rv_I_new[0], rv_I_new[1], rv_I_new[2], rv_I_new[3], rv_I_new[4], rv_I_new[5]);
    return 0; 
}