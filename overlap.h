#ifndef OVERLAP_H
#define OVERLAP_H

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <memory>

const double molecule_lidar_ratio = 8.0 * 3.1415926 / 3.0;

double get_diff(const std::unique_ptr<double[]> &ptr, const std::unique_ptr<double[]> &hptr, int n1, const int n2);
double runge_kutta(const std::unique_ptr<double[]> &signal_ratio, const std::unique_ptr<double[]> &molecule, const double, 
const std::unique_ptr<double[]> &distanceconst, int n, const double backscatter, const double, const double);
void openfile(std::ifstream &sin, const std::string &filename);

class bad_file{
    private:
        std::string filename;
    public:
        bad_file(const std::string &str) : filename(str) {}
        void mesg() {std::cout<<"The "<<filename<<" file does not exist!"<<std::endl;}
};
#endif