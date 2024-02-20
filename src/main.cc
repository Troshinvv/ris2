#include <iostream>
#include <TROOT.h>
#include <TSystem.h>

#define CALC_INCLUDE_PATH "/home/mikhail/ris2/src"
#define QnTools_LIB_PATH "QnTools::Base"
#define QnTools_INCLUDE_PATH "/home/mikhail/QnTools/install/lib/cmake/QnTools/../../../include/QnTools"


int main(int n_args, char** args) {
  if( n_args < 2 )
    throw std::runtime_error( "No arguments provided" );
  gInterpreter->Load(QnTools_LIB_PATH);
  gInterpreter->AddIncludePath(CALC_INCLUDE_PATH);
  gInterpreter->AddIncludePath(QnTools_INCLUDE_PATH);
  gROOT->Macro(args[1]);
}
