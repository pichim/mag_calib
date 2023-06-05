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

#include <stdint.h>

#pragma once

typedef struct magBiasEstimatorRLS_s {
    float p0;
    float lambda_min;
    float lambda;
    float b[3];
    float P[3][3];
    float Gamma[3][3];
} magBiasEstimatorRLS_t;

typedef struct magBiasEstimatorNLO_s {
    float Ts;
    float k[2];
    float b[3];
    float x[3];
} magBiasEstimatorNLO_t;

void magBiasEstimatorRLSInit(magBiasEstimatorRLS_t *magBiasEstimatorRLS, const float lambda_min, const float p0);
void magBiasEstimatorRLSReset(magBiasEstimatorRLS_t *magBiasEstimatorRLS);
void magBiasEstimatorRLSApply(magBiasEstimatorRLS_t *magBiasEstimatorRLS, float *mag, float *dmag, float *gyro);
void magBiasEstimatorRLSSolveRecursively(magBiasEstimatorRLS_t *magBiasEstimatorRLS, float *mag, float *dmag, float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign);

void magBiasEstimatorNLOInit(magBiasEstimatorNLO_t *magBiasEstimatorRLS, const float k1, const float k2, const uint32_t looptimeUs);
void magBiasEstimatorNLOReset(magBiasEstimatorNLO_t *magBiasEstimatorRLS);
void magBiasEstimatorNLOApply(magBiasEstimatorNLO_t *magBiasEstimatorRLS, float *mag, float *dmag, float *gyro);