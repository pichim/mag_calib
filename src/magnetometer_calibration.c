#include "magnetometer_calibration.h"

#include <stdint.h>

void magBiasEstimatorEstimatorInit(magBiasEstimatorEstimator_t *m)
{
    m->lambda_min = 0.9f;
    m->lambda = m->lambda_min;
    
    for (uint8_t i = 0; i < 3; i++) {
        m->b[i] = 0.0f;
        for (uint8_t j = 0; j < 3; j++) {
            if (i == j) {
                m->P[i][j] = 1.0e1f;
            } else {
                m->P[i][j] = 0.0f;
            }
        }
    }
}

void magBiasEstimatorEstimatorApply(magBiasEstimatorEstimator_t *m, float *mag, float *dmag, float *gyro)
{
    static float e[3];
    static float zn[3];
    static float Gamma[3][3];
    static float b[3];
    static float dP[3];

    // estimation error
    e[0] = dmag[0] + m->b[1] * gyro[2] - m->b[2] * gyro[1] + gyro[1] * mag[2] - gyro[2] * mag[1];
    e[1] = dmag[1] - m->b[0] * gyro[2] + m->b[2] * gyro[0] - gyro[0] * mag[2] + gyro[2] * mag[0];
    e[2] = dmag[2] + m->b[0] * gyro[1] - m->b[1] * gyro[0] + gyro[0] * mag[1] - gyro[1] * mag[0];
 
    // iteration 1
    zn[0] = 1.0f / ( m->lambda + gyro[2] * ( m->P[1][1] * gyro[2] - m->P[2][1] * gyro[1] ) - gyro[1] * (m->P[1][2] * gyro[2] - m->P[2][2] * gyro[1]) );
 
    Gamma[0][0] = zn[0] * ( m->P[0][2] * gyro[1] - m->P[0][1] * gyro[2] );
    Gamma[1][0] = zn[0] * ( m->P[1][2] * gyro[1] - m->P[1][1] * gyro[2] );
    Gamma[2][0] = zn[0] * ( m->P[2][2] * gyro[1] - m->P[2][1] * gyro[2] );

    zn[0] *= m->lambda;
 
    for (uint8_t i = 0; i < 3; i++) b[i] = m->b[i];
    m->b[0] += Gamma[0][0] * ( dmag[0] + b[1] * gyro[2] - b[2] * gyro[1] + gyro[1] * mag[2] - gyro[2] * mag[1] );
    m->b[1] += Gamma[1][0] * ( dmag[0] + b[1] * gyro[2] - b[2] * gyro[1] + gyro[1] * mag[2] - gyro[2] * mag[1] );
    m->b[2] += Gamma[2][0] * ( dmag[0] + b[1] * gyro[2] - b[2] * gyro[1] + gyro[1] * mag[2] - gyro[2] * mag[1] );
 
    dP[0] = m->P[1][0] * gyro[2] - m->P[2][0] * gyro[1];
    dP[1] = m->P[1][1] * gyro[2] - m->P[2][1] * gyro[1];
    dP[2] = m->P[1][2] * gyro[2] - m->P[2][2] * gyro[1];
    m->P[0][0] += Gamma[0][0] * dP[0];
    m->P[1][0] += Gamma[1][0] * dP[0];
    m->P[2][0] += Gamma[2][0] * dP[0];
    m->P[0][1] += Gamma[0][0] * dP[1];
    m->P[2][1] += Gamma[2][0] * dP[1];
    m->P[1][1] += Gamma[1][0] * dP[1];
    m->P[0][2] += Gamma[0][0] * dP[2];
    m->P[1][2] += Gamma[1][0] * dP[2];
    m->P[2][2] += Gamma[2][0] * dP[2];

    // iteration 2
    zn[1] = 1.0f / ( m->lambda + gyro[2] * ( m->P[0][0] * gyro[2] - m->P[2][0] * gyro[0] ) - gyro[0] * ( m->P[0][2] * gyro[2] - m->P[2][2] * gyro[0]) );
 
    Gamma[0][1] = zn[1] * ( m->P[0][0] * gyro[2] - m->P[0][2] * gyro[0] );
    Gamma[1][1] = zn[1] * ( m->P[1][0] * gyro[2] - m->P[1][2] * gyro[0] );
    Gamma[2][1] = zn[1] * ( m->P[2][0] * gyro[2] - m->P[2][2] * gyro[0] );

    zn[1] *= m->lambda;

    for (uint8_t i = 0; i < 3; i++) b[i] = m->b[i];
    m->b[0] += Gamma[0][1] * ( dmag[1] - b[0] * gyro[2] + b[2]*gyro[0] - gyro[0] * mag[2] + gyro[2] * mag[0] );
    m->b[1] += Gamma[1][1] * ( dmag[1] - b[0] * gyro[2] + b[2]*gyro[0] - gyro[0] * mag[2] + gyro[2] * mag[0] );
    m->b[2] += Gamma[2][1] * ( dmag[1] - b[0] * gyro[2] + b[2]*gyro[0] - gyro[0] * mag[2] + gyro[2] * mag[0] );

    dP[0] = m->P[2][0] * gyro[0] - m->P[0][0] * gyro[2];
    dP[1] = m->P[2][1] * gyro[0] - m->P[0][1] * gyro[2];
    dP[2] = m->P[2][2] * gyro[0] - m->P[0][2] * gyro[2];
    m->P[0][0] += Gamma[0][1] * dP[0];
    m->P[1][0] += Gamma[1][1] * dP[0];
    m->P[2][0] += Gamma[2][1] * dP[0];
    m->P[0][1] += Gamma[0][1] * dP[1];
    m->P[1][1] += Gamma[1][1] * dP[1];
    m->P[2][1] += Gamma[2][1] * dP[1];
    m->P[0][2] += Gamma[0][1] * dP[2];
    m->P[1][2] += Gamma[1][1] * dP[2];
    m->P[2][2] += Gamma[2][1] * dP[2];

    // iteration 3
    zn[2] = 1.0f / ( m->lambda + gyro[1] * ( m->P[0][0] * gyro[1] - m->P[1][0] * gyro[0] ) - gyro[0] * ( m->P[0][1]*gyro[1] - m->P[1][1]*gyro[0] ) );
 
    Gamma[0][2] = zn[2] * ( m->P[0][1] * gyro[0] - m->P[0][0] * gyro[1] );
    Gamma[1][2] = zn[2] * ( m->P[1][1] * gyro[0] - m->P[1][0] * gyro[1] );
    Gamma[2][2] = zn[2] * ( m->P[2][1] * gyro[0] - m->P[2][0] * gyro[1] );

    zn[2] *= m->lambda;

    for (uint8_t i = 0; i < 3; i++) b[i] = m->b[i];
    m->b[0] += Gamma[0][2] * ( dmag[2] + b[0] * gyro[1] - b[1] * gyro[0] + gyro[0] * mag[1] - gyro[1] * mag[0] );
    m->b[1] += Gamma[1][2] * ( dmag[2] + b[0] * gyro[1] - b[1] * gyro[0] + gyro[0] * mag[1] - gyro[1] * mag[0] );
    m->b[2] += Gamma[2][2] * ( dmag[2] + b[0] * gyro[1] - b[1] * gyro[0] + gyro[0] * mag[1] - gyro[1] * mag[0] );

    dP[0] = m->P[0][0] * gyro[1] - m->P[1][0] * gyro[0];
    dP[1] = m->P[0][1] * gyro[1] - m->P[1][1] * gyro[0];
    dP[2] = m->P[0][2] * gyro[1] - m->P[1][2] * gyro[0];
    m->P[0][0] += Gamma[0][2] * dP[0];
    m->P[1][0] += Gamma[1][2] * dP[0];
    m->P[2][0] += Gamma[2][2] * dP[0];
    m->P[0][1] += Gamma[0][2] * dP[1];
    m->P[1][1] += Gamma[1][2] * dP[1];
    m->P[2][1] += Gamma[2][2] * dP[1];
    m->P[0][2] += Gamma[0][2] * dP[2];
    m->P[1][2] += Gamma[1][2] * dP[2];
    m->P[2][2] += Gamma[2][2] * dP[2];

    // scale P with 1/lambda
    for (uint8_t i = 0; i < 3; i++) {
        for (uint8_t j = 0; j < 3; j++) {
            m->P[i][j] /= m->lambda;
        }
    }

    // adjust lambda
    m->lambda = m->lambda_min - ( m->lambda_min - 1.0f ) * ( zn[0] * zn[0] + zn[1] * zn[1] + zn[2] * zn[2] ) / 3.0f;

}