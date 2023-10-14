#include <iostream>
#include <fstream>
#include <chrono>

#include "axis.h"
#include "magnetometer_calibration.h"

#define LAMBDA_MIN 0.95f // minimal adaptive forgetting factor, range: [0.90, 0.99], currently tuned for 200 Hz
                         // (TASK_COMPASS_RATE_HZ) and update rate of compassBiasEstimatorApply(), not the mag readout
                         // rate, so it might need to be adjusted TASK_COMPASS_RATE_HZ is changed
#define P0 1.0e2f        // value to initialize P(0) = diag([P0, P0, P0]), typically in range: (1, 1000)

// run in terminal: .\output\mag_calib.exe "input/20231014_apex5_mag_on_tpu_03.txt"

int main(int argc, char *argv[])
{
    typedef struct mag_s {
        float magADC[XYZ_AXIS_COUNT];
    } mag_t;

    static mag_t mag;

    float data[6]; // mag, gyro
    FILE *dataInputFilePtr = fopen(argv[1], "r");
    std::string argInStr = argv[1];
    std::ofstream dataOutputFile("output/" + argInStr.substr(6, argInStr.length()));

    compassBiasEstimator_t compassBiasEstimator;
    compassBiasEstimatorInit(&compassBiasEstimator, LAMBDA_MIN, P0);
    compassBiasEstimatorUpdate(&compassBiasEstimator, LAMBDA_MIN, P0);

    while( fscanf(dataInputFilePtr, "%f, %f, %f, %f, %f, %f", &data[0], &data[1], &data[2], &data[3], &data[4], &data[5]) != EOF ) {
        
        for (int axis = 0; axis < XYZ_AXIS_COUNT; axis++) {
            mag.magADC[axis] = data[axis];
        }   

        // run compassBiasEstimatorApply and measure time
        std::chrono::steady_clock::time_point time_begin_ns = std::chrono::steady_clock::now();        
        compassBiasEstimatorApply(&compassBiasEstimator, mag.magADC);
        int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_begin_ns - std::chrono::steady_clock::now()).count();

        // print results to terminal and file
        //std::cout      << compassBiasEstimator.b[0] << ", " << compassBiasEstimator.b[1] << ", "
        //               << compassBiasEstimator.b[2] << ", " << compassBiasEstimator.lambda << ", " << time_elapsed_ns << std::endl;

        // print results to file
        dataOutputFile << compassBiasEstimator.b[0] << ", " << compassBiasEstimator.b[1] << ", "
                       << compassBiasEstimator.b[2] << ", " << compassBiasEstimator.lambda << std::endl;
    }

    fclose(dataInputFilePtr);
    dataOutputFile.close();

    return 0;
}