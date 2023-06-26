#include <iostream>
#include <fstream>
#include <chrono>
#include <stdio.h>

#include "magnetometer_calibration.h"

int main(int argc, char *argv[])
{
    float data[6]; // mag, gyro
    FILE *dataInputFilePtr = fopen("input/data_20230626_bf_LOG087.txt", "r");  

    std::ofstream dataOutputFile("output/mag_calib.txt");

    magBiasEstimatorRLS_t magBiasEstimatorRLS;
    magBiasEstimatorRLSInit(&magBiasEstimatorRLS, 0.99f, 1.0e0f);
    magBiasEstimatorNLO_t magBiasEstimatorNLO;
    magBiasEstimatorNLOInit(&magBiasEstimatorNLO, 20.0f, 0.2f, 20000);
    static float mag[3], dmag[3], gyro[3];
    static float magPast[3] = {0.0f, 0.0f, 0.0f};
    static float Ts = 0.1f;

    static bool is_first = true;

    while( fscanf(dataInputFilePtr, "%f, %f, %f, %f, %f, %f",
                                    &data[0], &data[1], &data[2], &data[3], &data[4], &data[5]) != EOF ) {
        
        mag[0]  = data[0];  mag[1] = data[1];  mag[2] = data[2];
        gyro[0] = data[3]; gyro[1] = data[4]; gyro[2] = data[5];

        if (is_first) {
            for (uint8_t i = 0; i < 3; i++) {
                magPast[i] = mag[i];
            }
            is_first = false;
        }
        
        for (uint8_t i = 0; i < 3; i++) {
            dmag[i] = ( mag[i] - magPast[i] ) / Ts;
            magPast[i] = mag[i];
        }

        std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
        
        magBiasEstimatorRLSApply(&magBiasEstimatorRLS, mag, dmag, gyro);
        magBiasEstimatorNLOApply(&magBiasEstimatorNLO, mag, dmag, gyro);

        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_begin).count();

        std::cout      << magBiasEstimatorRLS.b[0] << ", " << magBiasEstimatorRLS.b[1] << ", "
                       << magBiasEstimatorRLS.b[2] << ", " << magBiasEstimatorRLS.lambda << ", " 
                       << magBiasEstimatorNLO.b[0] << ", " << magBiasEstimatorNLO.b[1] << ", "
                       << magBiasEstimatorNLO.b[2] << ", "<< time_elapsed_ns << std::endl;
        dataOutputFile << magBiasEstimatorRLS.b[0] << ", " << magBiasEstimatorRLS.b[1] << ", "
                       << magBiasEstimatorRLS.b[2] << ", " << magBiasEstimatorRLS.lambda << ", " 
                       << magBiasEstimatorNLO.b[0] << ", " << magBiasEstimatorNLO.b[1] << ", "
                       << magBiasEstimatorNLO.b[2] << ", "<< time_elapsed_ns << std::endl;

    }

    fclose(dataInputFilePtr);
    dataOutputFile.close();

    return 0;
}