#include "magnetometer_calibration.h"

#include <stdint.h>

/**
 * Source and nomenclature:
 *   https://de.wikipedia.org/wiki/RLS-Algorithmus
 * Paper for estimation problem: https://www.roboticsproceedings.org/rss09/p50.pdf
 *   Adaptive Estimation of Measurement Bias in Three-Dimensional Field
 *   Sensors with Angular-Rate Sensors: Theory and Comparative Experimental Evaluation
 * Idea for adaptive forgetting factor is from:
 *   https://link.springer.com/book/10.1007/978-3-642-83530-8
 * and
 *   Ein Beitrag zur on-line adaptiven Regelung elektromechanischer Antriebsregelstrecken, Diss. 1997
 * Explicit Matrix inversion is avoided using the source:
 *   https://www.wiley.com/en-ie/Optimal+State+Estimation:+Kalman,+H+Infinity,+and+Nonlinear+Approaches-p-9780471708582
 *  
 * Problem formulation:
 * 
 *  armgin l2(e)
 *  e = y - y_hat
 *  Sw = Skew( gyro_x, gyro_y, gyro_z )
 *  y_hat = Sw * b, this is cross(gyro, b)
 *  y = d/dt mag + Sw * mag
 * 
 * Recursive Least Squares Algorithm with adaptive forgetting factor:
 * 
 *  for j = 1:3
 *  zn(j) = 1 / ( Sw(j,:) * P * Sw(j,:).' + lambda );
 *  Gamma(:,j) = P * Sw(j,:).' * zn(j);
 *  b = b + Gamma(:,j) * e(j);
 *  P = ( P - Gamma(:,j) * Sw(j,:) * P );
 *  end
 *  zn(j) = zn(j) * lambda;
 *  P = P / lambda;
 *  lambda = lambda_min + (1 - lambda_min) * ( zn.' * zn ) / 3.0;
 * 
 * Explicit form is:
 * 
 *  e = [dmag(i,1) + gyro(i,3) * ( b(2) - mag(i,2) ) - gyro(i,2) * ( b(3) - mag(i,3) ); ...
 *       dmag(i,2) - gyro(i,3) * ( b(1) - mag(i,1) ) + gyro(i,1) * ( b(3) - mag(i,3) ); ...
 *       dmag(i,3) + gyro(i,2) * ( b(1) - mag(i,1) ) - gyro(i,1) * ( b(2) - mag(i,2) )];
 * 
 *  // iteration 1
 *  dP = [P(2,1)*gyro(i,3) - P(3,1)*gyro(i,2); ...
 *        P(2,2)*gyro(i,3) - P(3,2)*gyro(i,2); ...
 *        P(3,2)*gyro(i,3) - P(3,3)*gyro(i,2)];
 *  zn(1) = 1 / ( lambda + gyro(i,3)*dP(2) - gyro(i,2)*dP(3) );
 *  g(1) = zn(1) * dP(1);
 *  g(2) = zn(1) * dP(2);
 *  g(3) = zn(1) * dP(3);
 *  b(1) = b(1) - e(1) * g(1);
 *  b(2) = b(2) - e(1) * g(2);
 *  b(3) = b(3) - e(1) * g(3);
 *  P(1,1) = P(1,1) - g(1) * dP(1);
 *  P(2,1) = P(2,1) - g(2) * dP(1);
 *  P(2,2) = P(2,2) - g(2) * dP(2);
 *  P(3,1) = P(3,1) - g(3) * dP(1);
 *  P(3,2) = P(3,2) - g(3) * dP(2);
 *  P(3,3) = P(3,3) - g(3) * dP(3);
 *
 *  // iteration 2
 *  dP = [P(3,1)*gyro(i,1) - P(1,1)*gyro(i,3); ...
 *        P(3,2)*gyro(i,1) - P(2,1)*gyro(i,3); ...
 *        P(3,3)*gyro(i,1) - P(3,1)*gyro(i,3)];
 *  zn(2) = 1/(lambda - gyro(i,3)*dP(1) + gyro(i,1)*dP(3));
 *  g(1) = zn(2) * dP(1);
 *  g(2) = zn(2) * dP(2);
 *  g(3) = zn(2) * dP(3);
 *  b(1) = b(1) - e(2) * g(1);
 *  b(2) = b(2) - e(2) * g(2);
 *  b(3) = b(3) - e(2) * g(3);
 *  P(1,1) = P(1,1) - g(1) * dP(1);
 *  P(2,1) = P(2,1) - g(2) * dP(1);
 *  P(2,2) = P(2,2) - g(2) * dP(2);
 *  P(3,1) = P(3,1) - g(3) * dP(1);
 *  P(3,2) = P(3,2) - g(3) * dP(2);
 *  P(3,3) = P(3,3) - g(3) * dP(3);
 *
 *  // iteration 3
 *  dP = [P(1,1)*gyro(i,2) - P(2,1)*gyro(i,1); ...
 *        P(2,1)*gyro(i,2) - P(2,2)*gyro(i,1); ...
 *        P(3,1)*gyro(i,2) - P(3,2)*gyro(i,1)];
 *  zn(3) = 1/(lambda + gyro(i,2)*dP(1) - gyro(i,1)*dP(2));
 *  g(1) = zn(3) * dP(1);
 *  g(2) = zn(3) * dP(2);
 *  g(3) = zn(3) * dP(3);
 *  b(1) = b(1) - e(3) * g(1);
 *  b(2) = b(2) - e(3) * g(2);
 *  b(3) = b(3) - e(3) * g(3);
 *  P(1,1) = P(1,1) - g(1) * dP(1);
 *  P(2,1) = P(2,1) - g(2) * dP(1);
 *  P(2,2) = P(2,2) - g(2) * dP(2);
 *  P(3,1) = P(3,1) - g(3) * dP(1);
 *  P(3,2) = P(3,2) - g(3) * dP(2);
 *  P(3,3) = P(3,3) - g(3) * dP(3);

 *  zn = zn * lambda;
 *  P = P / lambda;
 *  lambda = lambda_min + (1 - lambda_min) * ( zn.' * zn ) / 3.0;
 */
void compassBiasEstimatorInit(compassBiasEstimator_t *cBE, const float lambda_min, const float p0)
{
    compassBiasEstimatorUpdate(cBE, lambda_min, p0);
    
    cBE->lambda = lambda_min;

    // don't know if it is necessary to initialize those as 0
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

    // iteration 1: k = 0; i = 1; j = 2; sign =  1.0f;
    compassBiasEstimatorSolveIterative(cBE, e, zn, gyro, 0, 1, 2,  1.0f);
    // iteration 2: k = 1; i = 0; j = 2; sign = -1.0f;
    compassBiasEstimatorSolveIterative(cBE, e, zn, gyro, 1, 0, 2, -1.0f);
    // iteration 3: k = 2; i = 0; j = 1; sign =  1.0f;
    compassBiasEstimatorSolveIterative(cBE, e, zn, gyro, 2, 0, 1,  1.0f);

    for (uint8_t i = 0; i < 3; i++) {
        zn[i] *= cBE->lambda;
        for (uint8_t j = 0; j < 3; j++) {
            cBE->P[i][j] /= cBE->lambda;
        }
    }

    cBE->lambda = cBE->lambda_min - ( cBE->lambda_min - 1.0f ) * ( zn[0] * zn[0] + zn[1] * zn[1] + zn[2] * zn[2] ) / 3.0f;

}

void compassBiasEstimatorSolveIterative(compassBiasEstimator_t *cBE, const float *e, float *zn, const float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign)
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