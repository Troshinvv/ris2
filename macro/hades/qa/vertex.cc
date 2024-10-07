#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void vertex(){
  using namespace ris2;

  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style2D.cc" );

  auto plot = Plot( {2000, 1100} );

  auto au123 = Wrap<Histogram>( 
    "Au+Au E_{kin}=1.23A GeV"s, 
    "/home/mikhail/ris2/macro/hades/files/AuAu@123GEN9_100k.root"s,
    "event_header/vtx_z"s
  );
  au123.SetStyle( Style().SetMarker(-1).SetColor( kBlack ) );

  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 0.5, 1.0 } )
    .SetXAxis(Axis().SetTitle("vtx z (mm)").SetLo(-70.0).SetHi(10.0))
    .SetYAxis(Axis().SetTitle("counts").SetLo(1.0).SetHi(600.0))
    .AddToPlot( au123 )
  ;
  plot.AddSubPlot( std::vector<double>{ 0.5, 0.0, 1.0, 1.0 } )
    .AddToPlot( "/home/mikhail/ris2/macro/hades/files/AuAu@123GEN9_100k.root"s, "event_header/vtx_x_vtx_y"s )
    .SetXAxis(Axis().SetTitle("vtx x (mm)").SetLo(-5.0).SetHi(5.0))
    .SetYAxis(Axis().SetTitle("vtx y (mm)").SetLo(-5.0).SetHi(5.0))
    .SetZAxis(Axis().SetTitle("").SetLog())
  ;
  plot.Print( "/home/mikhail/ris2/macro/hades/qa/pictures/hades_vertex.png" );
};