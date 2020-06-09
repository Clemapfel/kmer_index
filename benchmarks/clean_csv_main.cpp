//
// Created by clem on 6/4/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <iostream>
#include <benchmarks/cleanup_csv.hpp>

int main(int argc, char** argv)
{
    std::cout << "cleaning " << argv[1] << "\n";
    cleanup_csv(argv[1]);
    std::cout << "done.\n";

}
