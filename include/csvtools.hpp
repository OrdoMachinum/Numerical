#ifndef CSVTOOLS_H
#define CSVTOOLS_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace csvtools{
    /* Writes a vector out to outFnm file */
    bool vecToCsv(std::string outFnm, const std::vector<double> & vecOfDoub);
    /* Writes a multivector to the outFnm file as a csv with delimiter delim */
    bool vecToCsv(std::string outFnm, const std::vector<std::vector<double> >& vecVecDouble, char delim);
    /* Reads out a csv into a multidimensional vector */
    bool csvToVects(std::string inFnm, std::vector<std::vector<double> > & vecVecDouble, char delim);
    /* Splits a dlimitet string into double values */
    int splitString(const std::string& inStr, std::vector<double> & outDoubleVect, char delim);
}


#endif