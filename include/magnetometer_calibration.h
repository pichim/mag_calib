/*
 * This file is part of Cleanflight and Betaflight.
 *
 * Cleanflight and Betaflight are free software. You can redistribute
 * this software and/or modify this software under the terms of the
 * GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Cleanflight and Betaflight are distributed in the hope that they
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.
 *
 * If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Source and nomenclature: https://de.wikipedia.org/wiki/RLS-Algorithmus
 * Paper for estimation problem: https://www.roboticsproceedings.org/rss09/p50.pdf
 *  Adaptive Estimation of Measurement Bias in Three-Dimensional Field
 *  Sensors with Angular-Rate Sensors: Theory and Comparative Experimental Evaluation
 * Idea for adaptive forgetting factor is from: https://link.springer.com/book/10.1007/978-3-642-83530-8
 * and                                        : Ein Beitrag zur on-line adaptiven Regelung elektromechanischer Antriebsregelstrecken, Diss. 1997
 * Explicit Matrix inversion is avoided using the source: https://www.wiley.com/en-ie/Optimal+State+Estimation:+Kalman,+H+Infinity,+and+Nonlinear+Approaches-p-9780471708582
 *
 * Problem formulation:
 * y_hat = Sw * b, Sw = Skew( gyro_x, gyro_y, gyro_z )
 * y     = d/dt mag + Sw * mag
 * e = y - y_hat, armgin l2(e)
 *
 * Recursive Least Squares Algorithm with adaptive forgetting factor:
 * for j = 1:3
 *      zn(j) = 1 / ( Sw(j,:) * P * Sw(j,:).' + lambda );
 *      Gamma(:,j) = P * Sw(j,:).' * zn(j);
 *      b = b + Gamma(:,j) * e(j);
 *      P = ( P - Gamma(:,j) * Sw(j,:) * P );
 *  end
 *  zn(j) = zn(j) * lambda;
 *  P = P / lambda;
 *  lambda = lambda_min + (1 - lambda_min) * ( zn.' * zn ) / 3.0;
*/

#include <stdint.h>

#pragma once

typedef struct compassBiasEstimator_s {
    float lambda_min, lambda, p0;
    float b[3];
    float P[3][3];
} compassBiasEstimator_t;

void compassBiasEstimatorInit(compassBiasEstimator_t *compassBiasEstimator, const float lambda_min, const float p0);
void compassBiasEstimatorUpdate(compassBiasEstimator_t *compassBiasEstimator, const float lambda_min, const float p0);
void compassBiasEstimatorApply(compassBiasEstimator_t *compassBiasEstimator, float *mag, const float *dmag, const float *gyro);
void compassBiasEstimatorSolveRecursively(compassBiasEstimator_t *compassBiasEstimator, const float *e, float *zn, const float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign);