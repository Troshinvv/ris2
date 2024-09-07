#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void pT_y_acceptance(){
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style2D.cc" );
  using namespace ris2;
  gStyle->SetPadTopMargin(0.15);

  auto plot = Plot( {2000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.5, 0.34, 1.0 } )
    .AddText( Text().SetText("TOF-400").SetPosition({0.2, 0.9}).SetSize(0.05) )
    .SetXAxis(Axis().SetTitle("y_{cm}"))
    .SetYAxis(Axis().SetTitle("p_{T} (GeV/c)"))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLog())
    .AddToPlot( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/qa.recent_vf_2024_06_20.root", "h2_y_pT_2212_tof400" )
    ;
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 0.34, 0.5 } )
    .AddText( Text().SetText("TOF-700").SetPosition({0.2, 0.9}).SetSize(0.05) )
    .SetXAxis(Axis().SetTitle("y_{cm}"))
    .SetYAxis(Axis().SetTitle("p_{T} (GeV/c)"))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLog())
    .AddToPlot( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/qa.recent_vf_2024_06_20.root", "h2_y_pT_2212_tof700" )
    ;
  plot.AddSubPlot( std::vector<double>{ 0.34, 0.0, 1.0, 1.0 } )
    .AddText( Text().SetText("Combined").SetPosition({0.2, 0.9}).SetSize(0.05) )
    .SetXAxis(Axis().SetTitle("y_{cm}"))
    .SetYAxis(Axis().SetTitle("p_{T} (GeV/c)"))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLog())
    .AddToPlot( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/qa.recent_vf_2024_06_20.root", "h2_y_pT_2212" )
    ;
  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/qa/pictures/pT_y_acceptance.png" );
};