#include <iostream>
#include <fstream>
#include <chrono>

#include "magnetometer_calibration.h"

#define XYZ_AXIS_COUNT 3

#define LAMBDA_MIN 0.99f // minimal adaptive forgetting factor
#define P0 1.0e0f        // value to initialize P(0) = diag([P0, P0, P0])

// run:.\output\mag_calib.exe "input/data_20230626_bf_LOG087.txt"

int main(int argc, char *argv[])
{
    float data[6]; // mag, gyro
    FILE *dataInputFilePtr = fopen(argv[1], "r");

    std::string argInStr = argv[1];
    std::ofstream dataOutputFile("output/mag_calib_" + argInStr.substr(argInStr.find_last_of("LOG")-2, argInStr.find_last_of("txt")));

    compassBiasEstimator_t compassBiasEstimator;
    compassBiasEstimatorInit(&compassBiasEstimator, LAMBDA_MIN, P0);
    compassBiasEstimatorUpdate(&compassBiasEstimator, LAMBDA_MIN, P0);

    typedef struct mag_s {
        float magADC[XYZ_AXIS_COUNT];
    } mag_t;

    mag_t mag, magPrevious;

    typedef enum {
        X = 0,
        Y,
        Z
    } axis_e;

    static const float dTime = 0.1f;
    static bool isFirstMagMeasurement = true;

    while( fscanf(dataInputFilePtr, "%f, %f, %f, %f, %f, %f",
                                    &data[0], &data[1], &data[2], &data[3], &data[4], &data[5]) != EOF ) {
        
        mag.magADC[X] = data[0]; mag.magADC[Y] = data[1]; mag.magADC[Z] = data[2];
        if (isFirstMagMeasurement) {
            isFirstMagMeasurement = false;
            magPrevious = mag;
        }
        const float dmag[3] =  { (mag.magADC[X] - magPrevious.magADC[X]) / dTime,
                                 (mag.magADC[Y] - magPrevious.magADC[Y]) / dTime,
                                 (mag.magADC[Z] - magPrevious.magADC[Z]) / dTime };
        magPrevious = mag;
        const float gyroADCfRadians[3] =  { data[3], data[4], data[5] };

        std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
        
        compassBiasEstimatorApply(&compassBiasEstimator, mag.magADC, dmag, gyroADCfRadians);

        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_begin).count();

        std::cout      << compassBiasEstimator.b[0] << ", " << compassBiasEstimator.b[1] << ", "
                       << compassBiasEstimator.b[2] << ", " << compassBiasEstimator.lambda << ", " << time_elapsed_ns << std::endl;
        dataOutputFile << compassBiasEstimator.b[0] << ", " << compassBiasEstimator.b[1] << ", "
                       << compassBiasEstimator.b[2] << ", " << compassBiasEstimator.lambda << ", " <<  time_elapsed_ns << std::endl;
    }

    fclose(dataInputFilePtr);
    dataOutputFile.close();

    return 0;
}