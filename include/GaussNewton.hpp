#ifndef GAUSSNEWTON_H
#define GAUSSNEWTON_H

#include <iostream>
#include <vector>

typedef std::vector<std::vector < double> > matrix_t;


class GaussNewton1D {
    public:
        /* Constructor with argument of a function pointer */
        GaussNewton1D(double (*targetFncPtr) (double xi, std::vector<double> params));

        //Static  members
        /* Calculates the sum of squares of the elements of vector v */
        static double calcSqrSum(std::vector<double> const & v);
        /* Returns the determinant of the matrix */
        static double calcDeterminant(const matrix_t & matr) ;
        /* Returns the j column of the matrix */
        static std::vector<double> getColumnOfMatrix(unsigned int const j, matrix_t const & m);        
        /* Returns the ij indexed submatrix of the m matrix */
        static matrix_t getSubMatrix (unsigned int i, unsigned int j, matrix_t const & m);
        /* Returns the transponated matrix */
        static matrix_t getTransponated (matrix_t const & m);
        /* Returns the adjungated matrix m */
        static matrix_t getAdjungated (matrix_t const & m) ;
        /* Returns the product of a matrix and a scalar */
        static matrix_t productWithNumber(matrix_t const & m, double c);
        /* Returns the product of two matrix */
        static matrix_t product(matrix_t const & leftM, matrix_t const & rightM);
        /* Returns the inverted matrix */
        static matrix_t getInverse(matrix_t const & m);
        /* Makes every element of m below absLowest zero */
        static matrix_t truncate(matrix_t const & m, double absLowest);
        /* Returns the scalar product of two vectors */
        static double scalarProduct (std::vector<double> const & a, std::vector<double> const & b);
        /* Returns the product of a matrix and a column vector */
        static std::vector<double> productMatrVect(matrix_t const & leftM, std::vector<double> const & rightV);
        /* Return the sum of two vector */
        static std::vector<double> sumVects (std::vector<double> const & a, std::vector<double> const & b);
        /* Returns the scaled vector with c scalar  */
        static std::vector<double> scaleVects (std::vector<double> const & v, double c);
        /* Returns the sum of two matrices */
        static matrix_t sumMatrcies (matrix_t const & A, matrix_t const & B);
        /* Prints m matrix to a stream */
        static void printMatrixToStream(matrix_t const & m, std::ostream outStream);

        // non-static members
        /* Returns the approximated partial derivative of the function repressented by targetFncPtr along the paramter with the index of paramIdx */
        double partDerivative(double xi, unsigned int paramIdx);
        /* Calculates the Jacobian of the error function */
        void calcJacobianOfErrorFunction();
        /* Calculates the next iteration of the parameter set, represented by the param vector */
        void calcNextBeta();
        

        // Setter getter functions
        /* The partial derivative calculation is done on a dinamic scale */
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