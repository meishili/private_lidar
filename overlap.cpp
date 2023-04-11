#include "overlap.h"

double get_diff(const std::unique_ptr<double[]> &ptr, const std::unique_ptr<double[]> &hptr, int n1, const int n2){
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    double temp4 = 0.0;
    double diff;
    int n = n2 - n1 + 1;
    for(int i = n1; i <= n2; i++){
        temp1 += hptr[i];
        temp2 += ptr[i];
        temp3 += ptr[i] * hptr[i];
        temp4 += hptr[i] * hptr[i];
    }
    diff = (double(n) * temp3 - temp1 * temp2) / (double(n) * temp4 - temp1 * temp1);
    return diff;
}

double runge_kutta(const std::unique_ptr<double[]> &signal_ratio, const std::unique_ptr<double[]> &molecule, const double lidar_ratio,
const std::unique_ptr<double[]> &distance, const int n, const double backscatter, const double lamda_0, const double lamda_r){
    double temp1, temp2, temp3, temp4, temp5, temp6;
    double k1, k2, k3, k4;

    temp1 = get_diff(molecule, distance, n - 3, n + 3) / molecule[n];
    temp2 = (1.0 - lamda_0 / lamda_r) * backscatter * lidar_ratio;
    temp3 = (1.0 - pow(lamda_0, 4) / pow(lamda_r, 4)) * molecule[n] * molecule_lidar_ratio;
    temp4 = get_diff(signal_ratio, distance, n - 3, n + 3) / signal_ratio[n];
    temp5 = molecule[n] + backscatter;
    temp6 = get_diff(molecule, distance, n - 3, n + 3);
    k1 = ((temp1 + temp2 + temp3 - temp4) * temp5 - temp6) * (distance[n] - distance[n - 1]);

    temp1 = get_diff(molecule, distance, n - 3, n + 2) / (molecule[n] + molecule[n - 1]) * 2.0;
    temp2 = (1.0 - lamda_0 / lamda_r) * (backscatter - k1 / 2.0) * lidar_ratio;
    temp3 = (1.0 - pow(lamda_0, 4) / pow(lamda_r, 4)) * (molecule[n - 1] + molecule[n]) / 2.0 * molecule_lidar_ratio;
    temp4 = get_diff(signal_ratio, distance, n - 3, n + 2) / (signal_ratio[n] + signal_ratio[n - 1]) * 2.0;
    temp5 = (molecule[n - 1] + molecule[n]) / 2.0 + backscatter - k1 / 2.0;
    temp6 = get_diff(molecule, distance, n - 3, n + 2);
    k2 = ((temp1 + temp2 + temp3 - temp4) * temp5 - temp6) * (distance[n] - distance[n - 1]);

    temp1 = get_diff(molecule, distance, n - 3, n + 2) / (molecule[n] + molecule[n - 1]) * 2.0;
    temp2 = (1.0 - lamda_0 / lamda_r) * (backscatter - k1 / 2.0) * lidar_ratio;
    temp3 = (1.0 - pow(lamda_0, 4) / pow(lamda_r, 4)) * (molecule[n - 1] + molecule[n]) / 2.0 * molecule_lidar_ratio;
    temp4 = get_diff(signal_ratio, distance, n - 3, n + 2) / (signal_ratio[n] + signal_ratio[n - 1]) * 2.0;
    temp5 = (molecule[n - 1] + molecule[n]) / 2.0 + backscatter - k2 / 2.0;
    temp6 = get_diff(molecule, distance, n - 3, n + 2);
    k3 = ((temp1 + temp2 + temp3 - temp4) * temp5 - temp6) * (distance[n] - distance[n - 1]);

    temp1 = get_diff(molecule, distance, n - 4, n + 2) / molecule[n];
    temp2 = (1.0 - lamda_0 / lamda_r) * (backscatter - k3) * lidar_ratio;
    temp3 = (1.0 - pow(lamda_0, 4) / pow(lamda_r, 4)) * molecule[n - 1] * molecule_lidar_ratio;
    temp4 = get_diff(signal_ratio, distance, n - 4, n + 2) / signal_ratio[n];
    temp5 = molecule[n - 1] + backscatter - k3;
    temp6 = get_diff(molecule, distance, n - 4, n + 2);
    k4 = ((temp1 + temp2 + temp3 - temp4) * temp5 - temp6) * (distance[n] - distance[n - 1]);

    double runge;
    runge = backscatter - (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    return runge;
}