#include <iostream>
#include <vector>
#include <chrono>
#include "csvtools.hpp"
#include "GaussNewton.hpp"

int main (int argc, char** argv) {
    
    
    matrix_t my;
    std::string fileName = "five.csv";

    if(argc == 2) {
        fileName = argv[1];
    }
    std::cout << "argc: " << argc << std::endl;


    if (!csvtools::csvToVects(fileName, ';',my)) {
        std::cout << "ERROR" << std::endl;
        return 1;
    }


    //auto out =  GaussNewton1D::getSubMatrix(1,1,my);
    //auto out =  GaussNewton1D::getTransponated(my);
    //std::cout << GaussNewton1D::calcDeterminant(out) << std::endl;
    auto start = std::chrono::steady_clock::now();
    //std::cout << GaussNewton1D::calcDeterminant(my) << std::endl;
    auto out = GaussNewton1D::truncate(GaussNewton1D::product(my, GaussNewton1D::getInverse(my)), 1E-14);
    //auto out = GaussNewton1D::product(my, my);
    //auto out = GaussNewton1D::getAdjungated(my);
    auto end = std::chrono::steady_clock::now();

    auto elapsedNanoSecs = end - start;

    std::cout << "Elapsed Time: " << elapsedNanoSecs.count() / 1.0E9 << " s " << std::endl;

    csvtools::vecToCsv("out.csv", ';', out);

    return 0;
}
