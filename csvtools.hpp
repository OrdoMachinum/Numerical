#ifndef CSVTOOLS_H
#define CSVTOOLS_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace csvtools{

    bool vecToCsv(std::string outFnm, const std::vector<double>& vecOfDoub);
    bool vecToCsv(std::string outFnm, char delim, const std::vector<std::vector<double> >& vecVecDouble);

    bool csvToVects(std::string inFnm, char delim, std::vector<std::vector<double> > & vecVecDouble);

    int splitString(const std::string& inStr, char delim, std::vector<double> & outDoubleVect);


}


#endif