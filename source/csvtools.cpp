#include "../include/csvtools.hpp"


bool csvtools::vecToCsv(std::string outFnm, const std::vector<double> & vecOfDoub) {
    std::ofstream outf(outFnm);
    if (!outf) {
        return false;
    }
    for (size_t i = 0; i < vecOfDoub.size(); i++)
    {
        outf << vecOfDoub.at(i) << std::endl;
    }
    outf.close();
    return true;
}

bool csvtools::vecToCsv(std::string outFnm, const std::vector<std::vector<double> > & vecOfDoub, char delim) {
    std::ofstream outf(outFnm);
    if (!outf) {
        return false;
    }
    //char delim = ';';
    std::string placeholder = "NaN";
    int numColumns = vecOfDoub.size();
    int maxLength = 0;
    std::vector<uint32_t> lengths;
    

    for (auto vec : vecOfDoub) {
       if (vec.size() > maxLength) {
           maxLength = vec.size();
       }
    }
    
    for (auto row :vecOfDoub) {
        for (auto c : row) {
            outf << c << delim;
        }
        outf << std::endl;
    }

    //std::cout << maxLength;

    outf.close();
    return true;

}

bool csvtools::csvToVects(std::string inFnm, char delim, std::vector<std::vector<double> > & vecVecDouble){
    std::ifstream in (inFnm);
    if(!in) {
        std::cout << "File reading error of " << inFnm << std::endl;
        return false;
    }
    std::string lineInput ="";
    std::vector<double> row;
    bool firstLine = true;
    while(!in.eof()) {
        std::getline(in,lineInput);
        if (lineInput.size() == 0) {
            continue;
        }
        splitString(lineInput, delim, row);
        vecVecDouble.push_back(row);
    }
    in.close();
    return true;
}

int csvtools::splitString(const std::string& inStr, char delim, std::vector<double> & outDoubleVect) {
     std::string word = "";
     outDoubleVect.clear();
     uint64_t startIndex = 0;
     uint64_t index = 0;
     while (index < inStr.size()) {
        if ((inStr[index] == delim) || (index >= inStr.size() - 1 )) {
            if(!word.empty()) {
                outDoubleVect.push_back(std::stod(word));
            }
            //std::cout << std::stod(word) <<std::endl;
            word.clear();
            ++index;
        } else {
            word.push_back(inStr[index]);
            ++index;
        }
     }
     return outDoubleVect.size();
}

