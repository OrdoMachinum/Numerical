#include "../include/GaussNewton.hpp"


GaussNewton1D::GaussNewton1D(double (*targetFncPtr) (double xi, std::vector<double> params)){
    m_targetFncPtr = targetFncPtr;
}

double GaussNewton1D::partDerivative(double xi, unsigned int paramIdx) { // numerical approx of partial derivative
    if (m_targetFncPtr == nullptr) {
        std::cout << " No Valid Function is loaded " << std::endl;
        return 0;
    }
    m_paramsVolatile = m_params;
    const double delta = m_params.at(paramIdx) * m_dinamicScale;

    m_paramsVolatile.at(paramIdx) = m_params.at(paramIdx) + delta;
    const double upper = m_targetFncPtr(xi, m_paramsVolatile);

    m_paramsVolatile.at(paramIdx) = m_params.at(paramIdx) - delta;
    const double lower = m_targetFncPtr(xi, m_paramsVolatile);

    return (upper-lower) / (2.0*delta);
    
}

void GaussNewton1D::calcJacobianOfErrorFunction() {
    if ((m_x.size() == 0) || (m_params.size() == 0) ) {
        std::cout << "Empty vectors" << std::endl;
        return;
    }
    m_jacobian.clear();
    m_jacobian.resize(m_x.size());
    for(unsigned int i = 0; i < m_x.size(); ++i) {
        for (unsigned int j = 0; j < m_params.size(); ++j) {
            m_jacobian.at(i).push_back( -1.0 * partDerivative(m_x.at(i), j)); // Calculation of Jacobian on an error function that is why -1 Partial derivtive
        }
    }
}

double GaussNewton1D::calcSqrSum(std::vector<double> const & v) {
    double ss {0};
    for(auto num : v) {
        ss += num*num;
    }
    return ss;
}

void GaussNewton1D::calcNextBeta() {
    calcJacobianOfErrorFunction();
    calcErrorVector();

    if ( (m_jacobian.size() == 0) || (m_r.size() == 0) || (m_params.size() == 0) ) {
        std::cout << "Not initalized members" << std::endl;
        return;
    }
    const auto prevParams = m_params;
    const auto transJac = getTransponated(m_jacobian);
    const auto inverted = getInverse(product(transJac, m_jacobian));

    const auto jacProdPred = product(inverted , transJac);

    m_params = sumVects(prevParams, scaleVects( productMatrVect(jacProdPred, m_r) , -1.0) );
}

double GaussNewton1D::calcDeterminant(const matrix_t & matr) {
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
                auto sm = GaussNewton1D::getSubMatrix(0, j, matr);
                d += sign * matr[0][j] * calcDeterminant(sm);
            }
    }
    return d;
    
}

matrix_t GaussNewton1D::getSubMatrix (unsigned int rowI, unsigned int columnJ, matrix_t const & m) {
    matrix_t subMatrix;
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

matrix_t GaussNewton1D::getTransponated(matrix_t const & m) {
    if(m.size() == 0) {
        std::cout << "Error, zero sized matrix cannot be transposed!! " <<std::endl;
        return m;
    }
    matrix_t mTransponated;
    std::vector<double> rowT;

    for(unsigned int j = 0; j < m.at(0).size(); ++j) {
        rowT.clear();
        for (unsigned int i = 0; i < m.size(); ++i) {
            rowT.push_back(m[i][j]);
        }
        mTransponated.push_back(rowT);
    }
    return mTransponated;
}

matrix_t GaussNewton1D::getAdjungated(matrix_t const & m) {
    if (m.size() == 0 ) {
        std::cout << "Error, zero sized matrix cannot be adjungated !! " <<std::endl;
        return m;
    }
    const matrix_t transposed = getTransponated(m);
    matrix_t adjungated;
    std::vector<double> rowAdj;
    short sign = 1;
    for(unsigned int i = 0; i < m.size(); ++i) {
        rowAdj.clear();
        sign = i % 2 == 0 ? -1 : 1;
        for(unsigned int j = 0; j < m.at(i).size(); ++j) {
            sign*=(-1);
            rowAdj.push_back(((double)sign)*calcDeterminant(getSubMatrix(i, j, transposed)));
        }
        adjungated.push_back(rowAdj);
    }
    return adjungated;
}
matrix_t GaussNewton1D::truncate(matrix_t const & m, double absLowest) {
    matrix_t truncated (m.size(), std::vector<double> (m.at(0).size()));
    for (unsigned int i = 0; i < m.size(); ++i) {
        for (unsigned int j = 0; j < m.at(i).size() ; ++j) {
            truncated[i][j] = (std::abs(m[i][j]) < absLowest ? 0.0 : m[i][j] );
        }
    }
    return truncated;
}

matrix_t GaussNewton1D::getInverse(matrix_t const & m) {
    return productWithNumber(getAdjungated(m), 1.0 / calcDeterminant(m));
}

matrix_t GaussNewton1D::productWithNumber(matrix_t const & m, double const c) {
    matrix_t scaled (m.size(), std::vector<double> (m[0].size())) ;
    for (unsigned int i = 0; i < m.size(); ++i) {
        for (unsigned int j = 0; j < m[0].size(); ++j) {
            scaled[i][j] = c * m[i][j];
        }
    }
    return scaled;
}

matrix_t GaussNewton1D::product(matrix_t const & leftM, matrix_t const & rightM) {
    if(leftM.at(0).size() != rightM.size() ){
        std::cout << " Matrix, matrix size mismatch for product " << std::endl;
        return leftM;
    }
    
    matrix_t prod (leftM.size(), std::vector<double>(rightM[0].size())) ;
    for (unsigned int i = 0; i < leftM.size(); ++i) {
        for (unsigned int j = 0; j < rightM[0].size(); ++j ) {
            prod[i][j] = scalarProduct(leftM[i], getColumnOfMatrix(j,rightM));
        }
    }
    return prod;
}

std::vector<double> GaussNewton1D::productMatrVect(matrix_t const & leftM, std::vector<double> const & rightV) {
    if(leftM.at(0).size() != rightV.size() ){
        std::cout << " Matrix, Vector size mismatch for product " << std::endl;
        return rightV;
    }
    std::vector<double> prod {};
    for (unsigned int i = 0; i < leftM.size(); ++i) {
        prod.push_back(scalarProduct(leftM[i], rightV));
    }
    return prod;
}
std::vector<double> GaussNewton1D::scaleVects (std::vector<double> const & v, double const c){
    std::vector<double> scaledV{};
    for(auto cell : v) {
        scaledV.push_back(cell*c);
    }
    return scaledV;
}

std::vector<double> GaussNewton1D::sumVects (std::vector<double> const & a, std::vector<double> const & b){
    if(a.size() != b.size()) {
        std::cout << "vector sum: size mismatch" << std::endl;
        return a;
    }
    std::vector<double> sum {};
    for(unsigned int i = 0; i < a.size(); ++i) {
        sum.push_back(a[i]+b[i]);
    }
    return sum;
}


double GaussNewton1D::scalarProduct (std::vector<double> const & a, std::vector<double> const & b) {
    if ( a.size() != b.size() ) {
        std::cout  << "Vector product: size mismatch" << std::endl;
        return 0;
    }
    double s = 0;
    for (unsigned int i = 0; i < a.size(); ++i) {
        s += a[i]*b[i];
    }
    return s;
}

std::vector<double> GaussNewton1D::getColumnOfMatrix(unsigned int const jColumn, matrix_t const & m) {
    std::vector<double> c {};
    if(jColumn >= m.at(0).size() ) {
        std::cout << "There is no column " << jColumn << " in this matrix." << std::endl;
        return c;
    }
    for (unsigned int i = 0; i < m.size(); ++i ) {
        c.push_back(m[i][jColumn]);
    }
    return c;
}


void GaussNewton1D::calcErrorVector() {
    if((m_x.size() == 0) || (m_y.size() == 0) || (m_params.size() == 0)) {
        std::cout << "Empty vector problem" << std::endl;
        return;
    }
    
    m_r.clear();
    for(unsigned int i =0 ; i < m_x.size(); ++i) {
        m_r.push_back(m_y.at(i) - m_targetFncPtr(m_x.at(i), m_params));
    }
}

void GaussNewton1D::setDinamicScale(double dinamicScale) {
     m_dinamicScale = dinamicScale;
}

double GaussNewton1D::getDinamicScale() const {
    return m_dinamicScale ;
}

void GaussNewton1D::setParams(std::vector<double> params){
    m_params = params;
}
std::vector<double> GaussNewton1D::getParams() const{
    return m_params;
}
void GaussNewton1D::setX(std::vector <double> x){
    m_x = x;
}
std::vector <double> GaussNewton1D::getX() const{
    return m_x;
}
void GaussNewton1D::setY(std::vector<double> raw){
    m_y = raw;
}
std::vector<double> GaussNewton1D::getY() const{
    return m_y;
}

std::vector<double> GaussNewton1D::getErrors() const{
    return m_r;
}

matrix_t GaussNewton1D::getJacbian() const {
    return m_jacobian;
}