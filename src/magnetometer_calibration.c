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
    static float dP[3];

    // estimation error
    e[0] = dmag[0] + m->b[1] * gyro[2] - m->b[2] * gyro[1] + gyro[1] * mag[2] - gyro[2] * mag[1];
    e[1] = dmag[1] - m->b[0] * gyro[2] + m->b[2] * gyro[0] - gyro[0] * mag[2] + gyro[2] * mag[0];
    e[2] = dmag[2] + m->b[0] * gyro[1] - m->b[1] * gyro[0] + gyro[0] * mag[1] - gyro[1] * mag[0];
 
    // iteration 1
    // i = 1
    // j = 2
    // k = 0
    // sign = 1
  //zn[k] = 1.0f / ( m->lambda + gyro[j] * ( m->P[i][i] * gyro[j] - m->P[j][i] * gyro[i] ) - gyro[i] * ( m->P[i][j] * gyro[j] - m->P[j][j] * gyro[i]) );
    zn[0] = 1.0f / ( m->lambda + gyro[2] * ( m->P[1][1] * gyro[2] - m->P[2][1] * gyro[1] ) - gyro[1] * ( m->P[1][2] * gyro[2] - m->P[2][2] * gyro[1]) );
 
  //Gamma[0][k] = sign * zn[k] * ( m->P[0][j] * gyro[i] - m->P[0][i] * gyro[j] );
  //Gamma[1][k] = sign * zn[k] * ( m->P[1][j] * gyro[i] - m->P[1][i] * gyro[j] );
  //Gamma[2][k] = sign * zn[k] * ( m->P[2][j] * gyro[i] - m->P[2][i] * gyro[j] );
    Gamma[0][0] = zn[0] * ( m->P[0][2] * gyro[1] - m->P[0][1] * gyro[2] );
    Gamma[1][0] = zn[0] * ( m->P[1][2] * gyro[1] - m->P[1][1] * gyro[2] );
    Gamma[2][0] = zn[0] * ( m->P[2][2] * gyro[1] - m->P[2][1] * gyro[2] );

  //zn[k] *= m->lambda;
    zn[0] *= m->lambda;
 
  //float db = dmag[k] + sign * ( b[i] * gyro[j] - b[j] * gyro[i] + gyro[i] * mag[j] - gyro[j] * mag[i] );
  //m->b[0] += Gamma[0][k] * db;
  //m->b[1] += Gamma[1][k] * db;
  //m->b[2] += Gamma[2][k] * db;
    float db = dmag[0] + m->b[1] * gyro[2] - m->b[2] * gyro[1] + gyro[1] * mag[2] - gyro[2] * mag[1];
    m->b[0] += Gamma[0][0] * db;
    m->b[1] += Gamma[1][0] * db;
    m->b[2] += Gamma[2][0] * db;
 
  //dP[0] = sign * ( m->P[i][0] * gyro[j] - m->P[j][0] * gyro[i] );
  //dP[1] = sign * ( m->P[i][1] * gyro[j] - m->P[j][1] * gyro[i] );
  //dP[2] = sign * ( m->P[i][2] * gyro[j] - m->P[j][2] * gyro[i] );
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
    // i = 0
    // j = 2
    // k = 1
    // sign = -1
  //zn[k] = 1.0f / ( m->lambda + gyro[j] * ( m->P[i][i] * gyro[j] - m->P[j][i] * gyro[i] ) - gyro[i] * ( m->P[i][j] * gyro[j] - m->P[j][j] * gyro[i]) );
    zn[1] = 1.0f / ( m->lambda + gyro[2] * ( m->P[0][0] * gyro[2] - m->P[2][0] * gyro[0] ) - gyro[0] * ( m->P[0][2] * gyro[2] - m->P[2][2] * gyro[0]) );
 
  //Gamma[0][k] = sign * zn[k] * ( m->P[0][j] * gyro[i] - m->P[0][i] * gyro[j] );
  //Gamma[1][k] = sign * zn[k] * ( m->P[1][j] * gyro[i] - m->P[1][i] * gyro[j] );
  //Gamma[2][k] = sign * zn[k] * ( m->P[2][j] * gyro[i] - m->P[2][i] * gyro[j] );
    Gamma[0][1] = zn[1] * ( m->P[0][0] * gyro[2] - m->P[0][2] * gyro[0] );
    Gamma[1][1] = zn[1] * ( m->P[1][0] * gyro[2] - m->P[1][2] * gyro[0] );
    Gamma[2][1] = zn[1] * ( m->P[2][0] * gyro[2] - m->P[2][2] * gyro[0] );

  //zn[k] *= m->lambda;
    zn[1] *= m->lambda;

  //float db = dmag[k] + sign * ( b[i] * gyro[j] - b[j] * gyro[i] + gyro[i] * mag[j] - gyro[j] * mag[i] );
  //m->b[0] += Gamma[0][k] * db;
  //m->b[1] += Gamma[1][k] * db;
  //m->b[2] += Gamma[2][k] * db;
    db = dmag[1] - m->b[0] * gyro[2] + m->b[2] * gyro[0] - gyro[0] * mag[2] + gyro[2] * mag[0];
    m->b[0] += Gamma[0][1] * db;
    m->b[1] += Gamma[1][1] * db;
    m->b[2] += Gamma[2][1] * db;

  //dP[0] = sign * ( m->P[i][0] * gyro[j] - m->P[j][0] * gyro[i] );
  //dP[1] = sign * ( m->P[i][1] * gyro[j] - m->P[j][1] * gyro[i] );
  //dP[2] = sign * ( m->P[i][2] * gyro[j] - m->P[j][2] * gyro[i] );
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
    // i = 0
    // j = 1
    // k = 2
    // sign = 1
  //zn[k] = 1.0f / ( m->lambda + gyro[j] * ( m->P[i][i] * gyro[j] - m->P[j][i] * gyro[i] ) - gyro[i] * ( m->P[i][j] * gyro[j] - m->P[j][j] * gyro[i]) );
    zn[2] = 1.0f / ( m->lambda + gyro[1] * ( m->P[0][0] * gyro[1] - m->P[1][0] * gyro[0] ) - gyro[0] * ( m->P[0][1]*gyro[1] - m->P[1][1]*gyro[0] ) );
 
  //Gamma[0][k] = sign * zn[k] * ( m->P[0][j] * gyro[i] - m->P[0][i] * gyro[j] );
  //Gamma[1][k] = sign * zn[k] * ( m->P[1][j] * gyro[i] - m->P[1][i] * gyro[j] );
  //Gamma[2][k] = sign * zn[k] * ( m->P[2][j] * gyro[i] - m->P[2][i] * gyro[j] );
    Gamma[0][2] = zn[2] * ( m->P[0][1] * gyro[0] - m->P[0][0] * gyro[1] );
    Gamma[1][2] = zn[2] * ( m->P[1][1] * gyro[0] - m->P[1][0] * gyro[1] );
    Gamma[2][2] = zn[2] * ( m->P[2][1] * gyro[0] - m->P[2][0] * gyro[1] );

  //zn[k] *= m->lambda;
    zn[2] *= m->lambda;

  //float db = dmag[k] + sign * ( b[i] * gyro[j] - b[j] * gyro[i] + gyro[i] * mag[j] - gyro[j] * mag[i] );
  //m->b[0] += Gamma[0][k] * db;
  //m->b[1] += Gamma[1][k] * db;
  //m->b[2] += Gamma[2][k] * db;
    db = dmag[2] + m->b[0] * gyro[1] - m->b[1] * gyro[0] + gyro[0] * mag[1] - gyro[1] * mag[0];
    m->b[0] += Gamma[0][2] * db;
    m->b[1] += Gamma[1][2] * db;
    m->b[2] += Gamma[2][2] * db;

  //dP[0] = sign * ( m->P[i][0] * gyro[j] - m->P[j][0] * gyro[i] );
  //dP[1] = sign * ( m->P[i][1] * gyro[j] - m->P[j][1] * gyro[i] );
  //dP[2] = sign * ( m->P[i][2] * gyro[j] - m->P[j][2] * gyro[i] );
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