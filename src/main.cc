#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

#include "process_buffered.hh"
#include "read_param.hh"

std::string Dir;
std::string Output;

float z;
float ext_alpha,ext_beta;

int 
main (int argc, char *argv[])
{
  Read_Param(argv[1]);
  ext_alpha = atof(argv[3]);
  ext_beta = atof(argv[4]);


  if (argv[2] == nullptr) {

    std::cerr << "Error! File for paths not given!\n";
    exit(EXIT_FAILURE);
  }

  std::ifstream inFile{argv[2]};

  if (inFile.is_open() == false) {

    std::cerr << "Error! File for paths cannot be opened!\n";
    exit(EXIT_FAILURE);
  }

  while ((inFile >> Dir >> z >> Output).eof() != true) {

    LIM_Buffered();

    std::cout << "\x1b[2K"
              << "\x1b[0G"
              << "Done for redshift: " << z << std::flush;
  }

  inFile.close();

  std::cout << "\x1b[2K"
            << "\x1b[0G"
            << "Done for all.\n"
            << std::flush;
}
