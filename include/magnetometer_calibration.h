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

#pragma once

typedef struct compassBiasEstimator_s {
    float lambda_min, lambda, p0;
    float b[3];
    float P[3][3];
} compassBiasEstimator_t;

void compassBiasEstimatorInit(compassBiasEstimator_t *compassBiasEstimator, const float lambda_min, const float p0);
void compassBiasEstimatorUpdate(compassBiasEstimator_t *compassBiasEstimator, const float lambda_min, const float p0);
void compassBiasEstimatorApply(compassBiasEstimator_t *compassBiasEstimator, float *mag, const float *dmag, const float *gyro);
void compassBiasEstimatorSolveIterative(compassBiasEstimator_t *compassBiasEstimator, const float *e, float *zn, const float *gyro, const unsigned k, const unsigned i, const unsigned j, const float sign);