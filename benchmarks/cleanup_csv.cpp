#include "cleanup_csv.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

// cleanup resulting csv for later parsing
void cleanup_csv(std::string path)
{
    std::ifstream file_in;
    file_in.open(path);

    if (file_in.fail())
    {
        std::cerr << "Error opening file " + path + "\naborting...";
        return;
    }

    // get timestamp

    std::string first_line;
    std::getline(file_in, first_line);

    std::replace(first_line.begin(), first_line.end(), ' ', '_');
    std::replace(first_line.begin(), first_line.end(), ':', '-');

    std::string out_path{path.begin(), path.end() - 7}; // - raw.csv
    out_path += first_line;
    out_path += ".csv";

    std::ofstream file_out;
    file_out.open(out_path);

    if (file_out.fail())
    {
        std::cerr << "Error opening file " + out_path + "\naborting...";
        return;
    }

    // discard and reformat header
    std::string line;
    bool past_header = false;

    while (std::getline(file_in, line))
    {
        if (not past_header)
        {
            // csv header line
            if (line.substr(0, 4) == "name")
            {
                line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());
                file_out << line + "\n";
                past_header = true;
            }
                // skip gbenchmark header
            else
                continue;
        }
        else
            file_out << line + "\n";
    }

    file_in.close();
    file_out.close();
}