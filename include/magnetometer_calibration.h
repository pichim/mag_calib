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

typedef struct magBiasEstimator_s {
    float lambda_min;
    float lambda;
    float b[3];
    float P[3][3];
} magBiasEstimator_t;

void magBiasEstimatorInit(magBiasEstimator_t *magBiasEstimator);
void magBiasEstimatorApply(magBiasEstimator_t *magBiasEstimator, float *mag, float *dmag, float *gyro);
void magBiasEstimatorSolveRecursively(magBiasEstimator_t *magBiasEstimator, float *mag, float *dmag, float *gyro, const uint8_t k, const uint8_t i, const uint8_t j, const float sign);