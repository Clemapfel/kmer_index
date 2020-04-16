//
// Created by Clemens Cords on 3/20/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <cassert>
#include <fstream>
#include <iostream>

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

    std::string out_path{path.begin(), path.end() - 3};
    out_path += "_clean.csv";

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

    std::cout << "[DEBUG] csv written to " + out_path + "\n";

    file_in.close();
    file_out.close();
}

int main(int argc, const char** argv)
{
    assert(argc == 2);
    cleanup_csv(argv[1]);

    std::cout << "finished processing csv at " << argv[1] << "\n";

}

