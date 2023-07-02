#include "magnetometer_calibration.h"

#include <stdint.h>

void magBiasEstimatorInit(magBiasEstimator_t *mBE, const float lambda_min, const float p0)
{
    mBE->lambda_min = lambda_min;
    mBE->p0 = p0;
    magBiasEstimatorReset(mBE);
}

void magBiasEstimatorReset(magBiasEstimator_t *mBE)
{
    mBE->lambda = mBE->lambda_min;
    for (uint8_t i = 0; i < 3; i++) {
        mBE->b[i] = 0.0f;
        for (uint8_t j = 0; j < 3; j++) {
            if (i == j) {
                mBE->P[i][j] = mBE->p0;
            } else {
                mBE->P[i][j] = 0.0f;
            }
        }
    }
}

void magBiasEstimatorApply(magBiasEstimator_t *mBE, float *mag, float *dmag, float *gyro)
{
    const float e[3] = { dmag[0] + gyro[2] * ( mBE->b[1] - mag[1] ) - gyro[1] * ( mBE->b[2] - mag[2] ),
                         dmag[1] - gyro[2] * ( mBE->b[0] - mag[0] ) + gyro[0] * ( mBE->b[2] - mag[2] ),
                         dmag[2] + gyro[1] * ( mBE->b[0] - mag[0] ) - gyro[0] * ( mBE->b[1] - mag[1] ) };

    // iteration 1 : k = 0; i = 1; j = 2; sign =  1.0f;
    magBiasEstimatorSolveRecursively(mBE, e, gyro, 0, 1, 2,  1.0f);
    // iteration 2 : k = 1; i = 0; j = 2; sign = -1.0f;
    magBiasEstimatorSolveRecursively(mBE, e, gyro, 1, 0, 2, -1.0f);
    // iteration 3 : k = 2; i = 0; j = 1; sign =  1.0f;
    magBiasEstimatorSolveRecursively(mBE, e, gyro, 2, 0, 1,  1.0f);

    for (uint8_t i = 0; i < 3; i++) {
        mBE->zn[i] *= mBE->lambda;
        for (uint8_t j = 0; j < 3; j++) {
            mBE->P[i][j] /= mBE->lambda;
        }
    }

    mBE->lambda = mBE->lambda_min - ( mBE->lambda_min - 1.0f ) * ( mBE->zn[0] * mBE->zn[0] + mBE->zn[1] * mBE->zn[1] + mBE->zn[2] * mBE->zn[2] ) / 3.0f;

}

void magBiasEstimatorSolveRecursively(magBiasEstimator_t *mBE, const float *e, float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign)
{
    const float dP[3] = { sign * ( mBE->P[0][i] * gyro[j] - mBE->P[0][j] * gyro[i] ),
                          sign * ( mBE->P[1][i] * gyro[j] - mBE->P[1][j] * gyro[i] ),
                          sign * ( mBE->P[2][i] * gyro[j] - mBE->P[2][j] * gyro[i] ) };

    mBE->zn[k] = 1.0f / ( mBE->lambda + sign * ( dP[i] * gyro[j] - dP[j] * gyro[i] ) );

    const float g[3] = { mBE->zn[k] * dP[0],
                         mBE->zn[k] * dP[1], 
                         mBE->zn[k] * dP[2] };

    for (uint8_t l = 0; l < 3; l++) {
        mBE->b[l] -= e[k] * g[l];
    }

    mBE->P[0][0] -= g[0] * dP[0];
    mBE->P[1][0] -= g[1] * dP[0];
    mBE->P[1][1] -= g[1] * dP[1];
    mBE->P[2][0] -= g[2] * dP[0];
    mBE->P[2][1] -= g[2] * dP[1];
    mBE->P[2][2] -= g[2] * dP[2];
    mBE->P[0][1] = mBE->P[1][0];
    mBE->P[0][2] = mBE->P[2][0];
    mBE->P[1][2] = mBE->P[2][1];
}