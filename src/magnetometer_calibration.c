#include "magnetometer_calibration.h"

#include <stdint.h>

void compassBiasEstimatorInit(compassBiasEstimator_t *cBE, const float lambda_min, const float p0)
{
    compassBiasEstimatorUpdate(cBE, lambda_min, p0);

    cBE->lambda = lambda_min;
    for (uint8_t i = 0; i < 3; i++) {
        cBE->b[i] = 0.0f;
        for (uint8_t j = 0; j < 3; j++) {
            if (i != j) {
                cBE->P[i][j] = 0.0f;
            }
        }
    }  
}

void compassBiasEstimatorUpdate(compassBiasEstimator_t *cBE, const float lambda_min, const float p0)
{
    cBE->lambda_min = lambda_min;
    cBE->p0 = p0;

    // update diagonal entries for faster convergence
    for (uint8_t i = 0; i < 3; i++) {
        for (uint8_t j = 0; j < 3; j++) {
            if (i == j) {
                cBE->P[i][j] = cBE->p0;
            }
        }
    } 
}

void compassBiasEstimatorApply(compassBiasEstimator_t *cBE, float *mag, const float *dmag, const float *gyro)
{
    const float e[3] = { dmag[0] + gyro[2] * ( cBE->b[1] - mag[1] ) - gyro[1] * ( cBE->b[2] - mag[2] ),
                         dmag[1] - gyro[2] * ( cBE->b[0] - mag[0] ) + gyro[0] * ( cBE->b[2] - mag[2] ),
                         dmag[2] + gyro[1] * ( cBE->b[0] - mag[0] ) - gyro[0] * ( cBE->b[1] - mag[1] ) };

    float zn[3];

    // iteration 1 : k = 0; i = 1; j = 2; sign =  1.0f;
    compassBiasEstimatorSolveRecursively(cBE, e, zn, gyro, 0, 1, 2,  1.0f);
    // iteration 2 : k = 1; i = 0; j = 2; sign = -1.0f;
    compassBiasEstimatorSolveRecursively(cBE, e, zn, gyro, 1, 0, 2, -1.0f);
    // iteration 3 : k = 2; i = 0; j = 1; sign =  1.0f;
    compassBiasEstimatorSolveRecursively(cBE, e, zn, gyro, 2, 0, 1,  1.0f);

    for (uint8_t i = 0; i < 3; i++) {
        zn[i] *= cBE->lambda;
        for (uint8_t j = 0; j < 3; j++) {
            cBE->P[i][j] /= cBE->lambda;
        }
    }

    cBE->lambda = cBE->lambda_min - ( cBE->lambda_min - 1.0f ) * ( zn[0] * zn[0] + zn[1] * zn[1] + zn[2] * zn[2] ) / 3.0f;

}

void compassBiasEstimatorSolveRecursively(compassBiasEstimator_t *cBE, const float *e, float *zn, const float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign)
{
    const float dP[3] = { sign * ( cBE->P[0][i] * gyro[j] - cBE->P[0][j] * gyro[i] ),
                          sign * ( cBE->P[1][i] * gyro[j] - cBE->P[1][j] * gyro[i] ),
                          sign * ( cBE->P[2][i] * gyro[j] - cBE->P[2][j] * gyro[i] ) };

    zn[k] = 1.0f / ( cBE->lambda + sign * ( dP[i] * gyro[j] - dP[j] * gyro[i] ) );

    const float g[3] = { zn[k] * dP[0],
                         zn[k] * dP[1], 
                         zn[k] * dP[2] };

    for (uint8_t l = 0; l < 3; l++) {
        cBE->b[l] -= e[k] * g[l];
    }

    cBE->P[0][0] -= g[0] * dP[0];
    cBE->P[1][0] -= g[1] * dP[0];
    cBE->P[1][1] -= g[1] * dP[1];
    cBE->P[2][0] -= g[2] * dP[0];
    cBE->P[2][1] -= g[2] * dP[1];
    cBE->P[2][2] -= g[2] * dP[2];
    cBE->P[0][1] = cBE->P[1][0];
    cBE->P[0][2] = cBE->P[2][0];
    cBE->P[1][2] = cBE->P[2][1];
}