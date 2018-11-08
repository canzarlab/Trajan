#include <vector>
#include <string>
#include "read_csv.h"
using namespace std;

/*
* Parses through csv file line by line and returns the data
* in vector of vector of strings.
*/
std::vector<std::vector<std::string> > CSVReader::getStringData(int from_row, int to_row)
{
	std::ifstream file(fileName);

	std::vector<std::vector<std::string> > dataList;

	std::string line = "";
	// Iterate through each line and split the content using delimeter
    int count = 0;
	while (getline(file, line))
	{
        if (count >= to_row) {
            break;
        }
		std::vector<std::string> vec;
        if (count >= from_row){
            boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
            dataList.push_back(vec);
        }
        count ++;
	}
	// Close the File
	file.close();

	return dataList;
}

std::vector<std::vector<double> > CSVReader::getDoubleFromStringData(int from_row, int to_row) {
    vector<vector<double> > data;
    std::vector<std::vector<string> > dataList = getStringData();
    for (int i = from_row; i < std::min(to_row, (int) dataList.size()); ++i)
    {
        vector<double> temp;
        for (const auto & u: dataList.at(i))
        {
            temp.push_back(std::stod(u));
        }
        data.push_back(temp);
    }
    return data;
}

std::vector<std::vector<double> > CSVReader::getDoubleData(int from_row, int to_row)
{
    vector<vector<double> > data;
    ifstream inputFile(fileName);
    int l = 0; // 0
    int count =0;
 
    while (inputFile) {
        if (count >= to_row) break;
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> record;
            int removefist = 0;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    if (removefist >= 0)    // matrix
                    record.push_back(stof(line));
                    else
                        removefist++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << fileName << " line " << l
                         << endl;
                    e.what();
                }
            }
            if (count >= from_row)
                data.push_back(record);
            count ++;
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << fileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}
