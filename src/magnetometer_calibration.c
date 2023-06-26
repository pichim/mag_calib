#include "magnetometer_calibration.h"

#include <stdint.h>

void magBiasEstimatorRLSInit(magBiasEstimatorRLS_t *mBE, const float lambda_min, const float p0)
{
    mBE->p0 = p0;
    mBE->lambda_min = lambda_min;
    magBiasEstimatorRLSReset(mBE);
}

void magBiasEstimatorRLSReset(magBiasEstimatorRLS_t *mBE)
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
            mBE->Gamma[i][j] = 0.0f;
        }
    }
}

void magBiasEstimatorRLSApply(magBiasEstimatorRLS_t *mBE, float *mag, float *dmag, float *gyro)
{
    // iteration 1 : k = 0; i = 1; j = 2; sign =  1.0f;
    magBiasEstimatorRLSSolveRecursively(mBE, mag, dmag, gyro, 0, 1, 2,  1.0f);
    // iteration 2 : k = 1; i = 0; j = 2; sign = -1.0f;
    magBiasEstimatorRLSSolveRecursively(mBE, mag, dmag, gyro, 1, 0, 2, -1.0f);
    // iteration 3 : k = 2; i = 0; j = 1; sign =  1.0f;
    magBiasEstimatorRLSSolveRecursively(mBE, mag, dmag, gyro, 2, 0, 1,  1.0f);
}

void magBiasEstimatorRLSSolveRecursively(magBiasEstimatorRLS_t *mBE, float *mag, float *dmag, float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign)
{
    static float zn[3];

    zn[k] = 1.0f / ( mBE->lambda + gyro[j] * ( mBE->P[i][i] * gyro[j] - mBE->P[j][i] * gyro[i] ) - gyro[i] * ( mBE->P[i][j] * gyro[j] - mBE->P[j][j] * gyro[i]) );
 
    for (uint8_t l = 0; l < 3; l++) mBE->Gamma[l][k] = sign * zn[k] * ( mBE->P[l][j] * gyro[i] - mBE->P[l][i] * gyro[j] );
  
    zn[k] *= mBE->lambda;

    float db = dmag[k] + sign * ( mBE->b[i] * gyro[j] - mBE->b[j] * gyro[i] + gyro[i] * mag[j] - gyro[j] * mag[i] );
    for (uint8_t l = 0; l < 3; l++) mBE->b[l] += mBE->Gamma[l][k] * db;

    float dP[3] = { sign * ( mBE->P[i][0] * gyro[j] - mBE->P[j][0] * gyro[i] ),
                    sign * ( mBE->P[i][1] * gyro[j] - mBE->P[j][1] * gyro[i] ),
                    sign * ( mBE->P[i][2] * gyro[j] - mBE->P[j][2] * gyro[i] ) };
    for (uint8_t l = 0; l < 3; l++) {
        for (uint8_t m = 0; m < 3; m++) {
            mBE->P[m][l] += mBE->Gamma[m][k] * dP[l];
        }
    }
    
    if (k == 2) {

        for (uint8_t l = 0; l < 3; l++) {
            for (uint8_t m = 0; m < 3; m++) {
                mBE->P[l][m] /= mBE->lambda;
            }
        }

        mBE->lambda = mBE->lambda_min + ( 1.0f - mBE->lambda_min ) * ( zn[0] * zn[0] + zn[1] * zn[1] + zn[2] * zn[2] ) / 3.0f;
    }
}