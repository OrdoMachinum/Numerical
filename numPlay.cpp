#include <iostream>
#include <vector>
#include <chrono>
#include "csvtools.hpp"
#include "GaussNewton.hpp"

int main () {
    
    
    std::vector<std::vector < double> > my;

    csvtools::csvToVects("five.csv", ';',my);


    auto out =  GaussNewton1D::getSubMatrix(1,1,my);

    //std::cout << GaussNewton1D::calcDeterminant(out) << std::endl;
    auto start = std::chrono::steady_clock::now();
    std::cout << GaussNewton1D::calcDeterminant(my) << std::endl;
    auto end = std::chrono::steady_clock::now();

    auto elapsedNanoSecs = end - start;

    std::cout << "Elapsed Time: " << elapsedNanoSecs.count() /1.0E9 << " s " << std::endl;

    csvtools::vecToCsv("out.csv", ';', out);

    return 0;
}
