#ifndef GAUSSNEWTON_H
#define GAUSSNEWTON_H

#include <iostream>
#include <vector>

typedef std::vector<std::vector < double> > matrix_t;


class GaussNewton1D {
    public:
        GaussNewton1D(double (*targetFncPtr) (double xi, std::vector<double> params));

        //Static  members
        static double calcSqrSum(std::vector<double> const & v);

        static double calcDeterminant(const matrix_t & matr) ;
        static std::vector<double> getColumnOfMatrix(unsigned int const j, matrix_t const & m);        
        
        static matrix_t getSubMatrix (unsigned int i, unsigned int j, matrix_t const & m);
        static matrix_t getTransponated (matrix_t const & m);
        static matrix_t getAdjungated (matrix_t const & m) ;
        static matrix_t productWithNumber(matrix_t const & m, double c);
        static matrix_t product(matrix_t const & a, matrix_t const & b);
        static matrix_t getInverse(matrix_t const & m);
        static matrix_t truncate(matrix_t const & m, double absLowest);
        static double scalarProduct (std::vector<double> const & a, std::vector<double> const & b);
        static std::vector<double> productMatrVect(matrix_t const & leftM, std::vector<double> const & rightV);
        static std::vector<double> sumVects (std::vector<double> const & a, std::vector<double> const & b);
        static std::vector<double> scaleVects (std::vector<double> const & v, double c);
        static matrix_t sumMatrcies (matrix_t const & A, matrix_t const & B);
        static void printMatrixToStream(matrix_t const & m, std::ostream outStream);

        // non-static members
        double partDerivative(double xi, unsigned int paramIdx);
        void calcJacobianOfErrorFunction();

        void calcNextBeta();

        double ErrorSumSquare();


        // Setter getter functions
        void setDinamicScale(double dinamicScale);
        double getDinamicScale() const;

        void setParams(std::vector<double> params);
        std::vector<double> getParams() const;

        void setX(std::vector <double> x);
        std::vector <double>  getX() const;

        void setY(std::vector<double> raw);
        std::vector<double> getY() const;

        void calcErrorVector();
        std::vector<double> getErrors() const;

        matrix_t getJacbian() const;
        
        
    private:
        std::vector <double> m_x; // Independent Variable
        std::vector <double> m_y; // Dependent Variable
        std::vector <double> m_r; // Error
        std::vector <double> m_params; // vector of actual parameters
        std::vector <double> m_paramsVolatile; // parameter vector for partial derivative calcuation
        matrix_t m_jacobian;

        double (*m_targetFncPtr) (double xi, std::vector<double> params); // function representing the function to fit

        double m_dinamicScale {0.05};
};


#endif