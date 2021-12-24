#include "GaussNewton.hpp"


GaussNewton1D::GaussNewton1D(double (*targetFncPtr) (double xi, std::vector<double> params)){
    _targetFncPtr = targetFncPtr;
}

double GaussNewton1D::partDerivative(double xi, unsigned int paramIdx) { // numerical approx of partial derivative
    if (_targetFncPtr == nullptr) {
        std::cout << " No Valid Function is loaded " << std::endl;
        return 0;
    }
    _paramsVolatile = _params;
    const double delta = _params.at(paramIdx) * _dinamicScale;

    _paramsVolatile.at(paramIdx) = _params.at(paramIdx) + delta;
    const double upper = _targetFncPtr(xi, _paramsVolatile);

    _paramsVolatile.at(paramIdx) = _params.at(paramIdx) - delta;
    const double lower = _targetFncPtr(xi, _paramsVolatile);

    return (upper-lower) / (2.0*delta);
    
}

void GaussNewton1D::calcJacobianOfErrorFunction() {
    if ((_x.size() == 0) || (_params.size() == 0) ) {
        std::cout << "Empty vectors" << std::endl;
        return;
    }
    _jacobian.clear();
    _jacobian.resize(_x.size());
    for(unsigned int i = 0; i < _x.size(); ++i) {
        for (unsigned int j = 0; j < _params.size(); ++j) {
            _jacobian.at(i).push_back( -1.0 * partDerivative(_x.at(i), j)); // Calculation of Jacobian on an error function that is why -1 Partial derivtive
        }
    }
}

void GaussNewton1D::calcNextBeta() {
    if ( (_jacobian.size() == 0) || (_r.size() == 0) || (_params.size() == 0) ) {
        std::cout << "Not initalized members" << std::endl;
        return;
    }
    
}

double GaussNewton1D::calcDeterminant(const std::vector<std::vector<double> > & matr) {
    const unsigned int sizeDet { matr.size() };
    if (sizeDet == 0) {
        std::cout << " empty matrix " <<std::endl;
        return 0;
    }
    for(unsigned int i = 0; i < matr.size(); ++i) {
        if (matr.at(i).size() != matr.size()) {
            std::cout << " NON RECTANGULAR MATRIX !! " << std::endl;
            return -1;
        }
    }
    double d {0.0};
    
    switch (sizeDet) {
        case 1:
            return matr[0][0];
            break;
        case 2:
            return matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0];
            break;
        default:
            short sign = -1;
            for (unsigned int j = 0; j < sizeDet; ++j) {
                sign *= -1;
                auto sm = getSubMatrix(0, j, matr);
                d += sign * matr[0][j] * calcDeterminant(sm);
            }
    }
    return d;
    
}

std::vector<std::vector < double> > GaussNewton1D::getSubMatrix (unsigned int rowI, unsigned int columnJ, std::vector<std::vector<double> > const & m) {
    std::vector<std::vector < double> > subMatrix;
    for(unsigned int i = 0; i < m.size(); ++i) {
        if(i==rowI){
            continue;
        }
        std::vector<double> rowActual;
        for(unsigned int j = 0; j < m.at(i).size(); ++j) {
            if(j==columnJ) {
                continue;
            }
            rowActual.push_back(m[i][j]);
        }
        subMatrix.push_back(rowActual);
    }
    return subMatrix;
}




void GaussNewton1D::calcErrorVector() {
    if((_x.size() == 0) || (_y.size() == 0) ) {
        std::cout << "Empty vector problem" << std::endl;
        return;
    }
    _r.clear();
    for(unsigned int i =0 ; i < _x.size(); ++i) {
        _r.push_back(_y.at(i) - _targetFncPtr(_x.at(i), _params));
    }
}

void GaussNewton1D::setDinamicScale(double dinamicScale) {
     _dinamicScale = dinamicScale;
}

double GaussNewton1D::getDinamicScale() const {
    return _dinamicScale ;
}