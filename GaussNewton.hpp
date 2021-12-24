#ifndef GAUSSNEWTON_H
#define GAUSSNEWTON_H

#include <iostream>
#include <vector>


class GaussNewton1D {
    public:
        GaussNewton1D(double (*targetFncPtr) (double xi, std::vector<double> params));

        double partDerivative(double xi, unsigned int paramIdx);
        void calcJacobianOfErrorFunction();
        static double calcDeterminant(const std::vector<std::vector<double> > & matr) ;

        static std::vector<std::vector < double> > getSubMatrix (unsigned int i, unsigned int j, std::vector<std::vector<double> > const & m);

        void calcNextBeta();

        double ErrorSumSquare();

        void setDinamicScale(double dinamicScale);
        double getDinamicScale() const;

        void setParams(std::vector<double> params);
        std::vector<double> getParams() const;

        void setX(std::vector <double> x);
        double getX() const;

        void setY(std::vector<double> params);
        std::vector<double> getY() const;

        void calcErrorVector();
        std::vector<double> getErrors() const;
        
        
    private:
        std::vector <double> _x; // Independent Variable
        std::vector <double> _y; // Dependent Variable
        std::vector <double> _r; // Error
        std::vector <double> _params; // vector of actual parameters
        std::vector <double> _paramsVolatile; // parameter vector for partial derivative calcuation
        std::vector<std::vector<double> > _jacobian;

        double (*_targetFncPtr) (double xi, std::vector<double> params);

        double _dinamicScale {0.05};

        

};


#endif