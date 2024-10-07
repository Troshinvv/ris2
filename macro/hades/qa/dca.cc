#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void dca(){
  using namespace ris2;

  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );

  auto plot = Plot( {1000, 1100} );

  auto au123 = Wrap<Histogram>( 
    "Au+Au E_{kin}=1.23A GeV"s, 
    "/home/mikhail/ris2/macro/hades/files/AuAu@123GEN9_100k.root"s,
    "mdc_vtx_tracks/dca_xy"s
  );
  au123.SetStyle( Style().SetMarker(-1).SetColor( kBlack ) );

  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(Axis().SetTitle("DCA (mm)").SetLo(-15.0).SetHi(15.0))
    .SetYAxis(Axis().SetTitle("counts").SetLo(1.0).SetHi(30000.0))
    .AddToPlot( au123 )
  ;
  plot.Print( "/home/mikhail/ris2/macro/hades/qa/pictures/hades_dca.png" );
};