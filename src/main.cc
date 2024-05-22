#include <iostream>
#include <TROOT.h>
#include <TSystem.h>

#define CALC_INCLUDE_PATH "/home/mikhail/ris2/src"

int main(int n_args, char** args) {
  if( n_args < 2 )
    throw std::runtime_error( "No arguments provided" );
  gInterpreter->AddIncludePath(CALC_INCLUDE_PATH);
  gROOT->Macro(args[1]);
}
