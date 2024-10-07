#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void histo2d_example(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  using namespace ris2;

  auto plot = Plot( {1000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.05).SetHi(1.1))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.5))
      .AddToPlot( "/home/mikhail/bmn_run8/efficiency.2024.04.03.root", "efficiency_2212_tof" )
      ;
  plot.Print( "/home/mikhail/ris2/macro/pictures/histo2d_example.png" );
};