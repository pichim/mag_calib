#include "magnetometer_calibration.h"

#include <string.h>

#include "maths.h"

// initialize the compass bias estimator
void compassBiasEstimatorInit(compassBiasEstimator_t *cBE, const float lambda_min, const float p0)
{
    memset(cBE, 0, sizeof(*cBE)); // zero contained IEEE754 floats
    // create identity matrix
    for (unsigned i = 0; i < 4; i++) {
        cBE->U[i][i] = 1.0f;
    }

    compassBiasEstimatorUpdate(cBE, lambda_min, p0); 

    cBE->lambda = lambda_min;
}

// reset / update the compass bias estimator, this can be used after the compass bias estimator did
// already run to achieve faster convergence for the next run
void compassBiasEstimatorUpdate(compassBiasEstimator_t *cBE, const float lambda_min, const float p0)
{
    cBE->lambda_min = lambda_min;
    // update diagonal entries for faster convergence
    for (unsigned i = 0; i < 4; i++) {
        cBE->D[i] = p0;
    } 
}

// apply one estimation step of the compass bias estimator
void compassBiasEstimatorApply(compassBiasEstimator_t *cBE, float *mag, const float scaleMag)
{
    float magScaled[3] = {mag[0] / scaleMag, mag[1] / scaleMag, mag[2] / scaleMag};

    // update phi
    float phi[4];
    phi[0] = sq(magScaled[0]) + sq(magScaled[1]) + sq(magScaled[2]);
    for (unsigned i = 0; i < 3; i++) {
        phi[i + 1] = magScaled[i];
    }

    // update e
    float e = 1.0f;
    for (unsigned i = 0; i < 4; i++) {
        e -= phi[i] * cBE->theta[i];
    }

    // U D U^T
    float f[4];
    float v[4];
    for (unsigned i = 0; i < 4; i++) {
        f[i] = 0.0f;
        for (unsigned j = 0; j <= i; j++) {
            f[i] += cBE->U[j][i] * phi[j];
        }
        v[i] = cBE->D[i] * f[i];
    }

    // first iteration
    float alpha[4];
    float k[4] = {0};
    alpha[0] = cBE->lambda + v[0] * f[0];
    cBE->D[0] /= alpha[0];
    k[0] = v[0];

    // rest of the iterations
    for (unsigned i = 1; i < 4; i++) {
        alpha[i] = alpha[i - 1] + v[i] * f[i];
        cBE->D[i] *= alpha[i - 1] / (alpha[i] * cBE->lambda);
        for (unsigned j = 0; j < i; j++) {
            float dU = -(f[i] / alpha[i - 1]) * k[j];
            k[j] += v[i] * cBE->U[j][i];
            cBE->U[j][i] += dU;
        }
        k[i] += v[i];
    }

    // parameter-update
    for (unsigned i = 0; i < 4; i++) {
        cBE->theta[i] += (k[i] / alpha[3]) * e;
    }

    // bias update
    for (unsigned i = 0; i < 3; i++) {
        cBE->b[i] = (-0.5f * cBE->theta[i + 1] / cBE->theta[0]) * scaleMag;
    }

    // compute zn
    float U_v;
    float phiTrans_U_v = 0.0f;
    for (unsigned i = 0; i < 4; i++) {
        U_v = 0.0f;
        for (unsigned j = i; j < 4; j++) {
            U_v += cBE->U[i][j] * v[j];
        }
        phiTrans_U_v += phi[i] * U_v;
    }
    float zn = cBE->lambda / (cBE->lambda + phiTrans_U_v);

    // update lambda
    cBE->lambda = cBE->lambda_min + (1.0f - cBE->lambda_min) * sq(zn);
}