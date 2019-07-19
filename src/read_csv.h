#ifndef READCSV_H
#define READCSV_H
/*
 * A class to read data from a csv file.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
class CSVReader
{
	std::string fileName;
	std::string delimeter;

public:
	CSVReader(std::string filename, std::string delm = ",") :
			fileName(filename), delimeter(delm)
	{ }

	// Function to fetch data from a CSV File
	std::vector<std::vector<std::string> > getStringData(int from_row = 0, int to_row = 1000000);
    
    std::vector<std::vector<double> > getDoubleFromStringData(int from_row = 0, int to_row = 1000000);
    
    std::vector<std::vector<double> > getDoubleData(int from_row = 0, int to_row = 10000);
};

#endif
