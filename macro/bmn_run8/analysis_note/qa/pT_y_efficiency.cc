#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void pT_y_efficiency(){
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style2D.cc" );
  using namespace ris2;

  auto plot = Plot( {1100, 1000} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(Axis().SetTitle("y_{cm}"))
    .SetYAxis(Axis().SetTitle("p_{T} (GeV/c)"))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLo(0.005).SetHi(1.0))
    .AddToPlot( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/efficiency.2024.04.03.root", "efficiency_2212_tof" )
    ;
  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/qa/pictures/pT_y_efficiency.png" );
};