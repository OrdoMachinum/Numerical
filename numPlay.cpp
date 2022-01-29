#include <iostream>
#include <vector>
#include <chrono>
#include "include/csvtools.hpp"
#include "include/GaussNewton.hpp"

double funcToFit (double x, std::vector<double> params) {
    return params[0] * (x-params[1])*(x-params[1]) + params[2];
}

int main (int argc, char** argv) {
    
    GaussNewton1D fitter(funcToFit);
    fitter.setParams({-0.01, 3, 300});

    matrix_t my;
    std::string fileName = "noised.csv";

    if(argc == 2) {
        fileName = argv[1];
    }
    std::cout << "argc: " << argc << std::endl;


    if (!csvtools::csvToVects(fileName, my, ';')) {
        std::cout << "ERROR" << std::endl;
        return 1;
    }

    fitter.setX(GaussNewton1D::getColumnOfMatrix(0, my));
    fitter.setY(GaussNewton1D::getColumnOfMatrix(1, my));
    
    matrix_t params{};
    fitter.setDinamicScale(0.02);

    //auto out =  GaussNewton1D::getSubMatrix(1,1,my);
    //auto out =  GaussNewton1D::getTransponated(my);
    //std::cout << GaussNewton1D::calcDeterminant(out) << std::endl;    
    auto start = std::chrono::steady_clock::now();
    //std::cout << GaussNewton1D::calcDeterminant(my) << std::endl;
    //auto out = GaussNewton1D::truncate(GaussNewton1D::product(my, GaussNewton1D::getInverse(my)), 1E-14);
    //auto out = GaussNewton1D::product(my, my);
    //auto out = GaussNewton1D::getAdjungated(my);
    fitter.calcErrorVector();
    fitter.calcJacobianOfErrorFunction();
    params.push_back(fitter.getParams());
    std::cout << "Err square sum: " << GaussNewton1D::calcSqrSum( fitter.getErrors() ) << std::endl;
    for(unsigned int i = 0 ; i < 30 ; ++i){
        fitter.calcNextBeta();
        params.push_back(fitter.getParams());
        std::cout << "Err square sum: " << GaussNewton1D::calcSqrSum( fitter.getErrors() ) << std::endl;
    }
    
    auto end = std::chrono::steady_clock::now();

    auto elapsedNanoSecs = end - start;

    std::cout << "Elapsed Time: " << elapsedNanoSecs.count() / 1.0E9 << " s " << std::endl;

    csvtools::vecToCsv("out.csv", fitter.getErrors());
    csvtools::vecToCsv("outJac.csv", fitter.getJacbian(),';');

    csvtools::vecToCsv("params.csv", params, '\t');
    
    std::cout << "END" << std::endl;
    return 0;
}
