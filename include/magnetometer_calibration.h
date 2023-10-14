#pragma once

typedef struct compassBiasEstimator_s {
    float lambda_min, lambda;
    float b[3];
    float theta[4];
    float U[4][4];
    float D[4];
} compassBiasEstimator_t;

void compassBiasEstimatorInit(compassBiasEstimator_t *compassBiasEstimator, const float lambda_min, const float p0);
void compassBiasEstimatorUpdate(compassBiasEstimator_t *compassBiasEstimator, const float lambda_min, const float p0);
void compassBiasEstimatorApply(compassBiasEstimator_t *cBE, float *mag);