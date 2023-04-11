#include <iostream>
#include "overlap.h"

int main(int argc, char *argv[]){
    int size = std::stoi(argv[1]);
    std::unique_ptr<double[]> distance(new double[size]);
    std::unique_ptr<double[]> p532(new double[size]);
    std::unique_ptr<double[]> p607(new double[size]);
    std::unique_ptr<double[]> signal_ratio(new double[size]);
    std::unique_ptr<double[]> molecule(new double[size]);
    std::ifstream mol_in;
    mol_in.open("air.txt");
    std::ifstream signal_in;
    signal_in.open("signal.txt");
    for(int i = 0; i < size; i++){
        mol_in>>molecule[i];
        molecule[i] = molecule[i] / molecule_lidar_ratio;
        signal_in>>distance[i]>>p532[i]>>p607[i];
        distance[i] /= 1000.0;
        signal_ratio[i] = p607[i] / p532[i];
    }
    mol_in.close();
    signal_in.close();
    double resolution = distance[1] - distance[0];
    double lidar_ratio;
    std::unique_ptr<double[]> raman_backscatter(new double[size]);
    std::unique_ptr<double[]> fernald_backscatter(new double[size]);
    int reference = size / 20 * 3 + 1;
    raman_backscatter[reference] = molecule[reference];
    fernald_backscatter[reference] = molecule[reference];
    std::cout<<fernald_backscatter[reference]<<std::endl;
    double lamda_0, lamda_r;
    std::cout<<"please enter wavelength of elastic and Raman ";
    std::cin>>lamda_0>>lamda_r;
    std::unique_ptr<std::unique_ptr<double[]>[]> overlap(new std::unique_ptr<double[]>[100]);
    for(int i = 0; i < 100; i++){
        overlap[i] = std::make_unique<double[]>(size);
    }
    std::unique_ptr<double[]> x = std::make_unique<double[]>(size);
    x[reference] = p532[reference] * distance[reference] * distance[reference];
    double error = 0.0;
    double error_sum = 100000000.0;
    double min_lidar_ratio;
    std::ofstream overlap_out;
    overlap_out.open("overlap.txt");
    std::ofstream backscatter_out;
    backscatter_out.open("backscatter.txt");
    for(int i = 0; i < 100; i++){
        lidar_ratio = double(i) + 1.0;
        double sum = 0.0;
        double temp1, temp2, temp3;
        for(int j = reference - 1; j > 0; j--){
            raman_backscatter[j] = runge_kutta(signal_ratio, molecule, lidar_ratio, distance, j + 1, raman_backscatter[j + 1], lamda_0, lamda_r);
            sum += (raman_backscatter[j] + raman_backscatter[j + 1]) * (1.0 + lamda_0 / lamda_r) * lidar_ratio;
            sum += (molecule[j] + molecule[j + 1]) * (1.0 + pow(lamda_0, 4) / pow(lamda_r, 4)) * molecule_lidar_ratio;
            overlap[i][j] = p607[j] * distance[j] * distance[j] * molecule[reference] / 
                p607[reference] / distance[reference] / distance[reference] / molecule[j];
            overlap[i][j] *= exp(-resolution * sum / 2.0);
            x[j] = p532[j] * distance[j] * distance[j] / overlap[i][j];
            temp1 = x[j] * exp((lidar_ratio - molecule_lidar_ratio) * (molecule[j] + molecule[j + 1]) * resolution);
            temp2 = x[j + 1] / (fernald_backscatter[j + 1] + molecule[j + 1]);
            temp3 = lidar_ratio * (x[j + 1] + temp1) * resolution;
            fernald_backscatter[j] = temp1 / (temp2 + temp3) - molecule[j];
            error += pow((raman_backscatter[j] - fernald_backscatter[j]), 2);
        }
        for(int j = 0; j < reference; j++){
            overlap_out<<lidar_ratio<<"  "<<distance[j]<<"  "<<overlap[i][j]<<std::endl;
            backscatter_out<<lidar_ratio<<"  "<<distance[j]<<"  "<<fernald_backscatter[j]<<"  "<<raman_backscatter[j]<<std::endl;
        }
        if(error < error_sum){
            min_lidar_ratio = lidar_ratio;
            error_sum = error;
        }
    }
    overlap_out.close();
    std::cout<<"lidar ratio "<<min_lidar_ratio<<std::endl;
    return 0;
}