#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void vtx_distribution(){
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style2D.cc" );
  using namespace ris2;

  auto plot = Plot( {3000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 0.33, 1.0 } )
    .SetXAxis(Axis().SetTitle("x_{vtx} (mm)").SetLo(-5.0).SetHi(5.0))
    .SetYAxis(Axis().SetTitle("N tracks").SetLo(0).SetHi(500))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLo(1).SetHi(280) )
    .AddToPlot( "/home/mikhail/bmn_run8/qa.full_vf_2024_07_08.root", "h2_vtx_x_multiplicity" )
    ;
  plot.AddSubPlot( std::vector<double>{ 0.33, 0.0, 0.66, 1.0 } )
    .SetXAxis(Axis().SetTitle("y_{vtx} (mm)").SetLo(-5.0).SetHi(5.0))
    .SetYAxis(Axis().SetTitle("N tracks").SetLo(0).SetHi(500))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLog())
    .AddToPlot( "/home/mikhail/bmn_run8/qa.full_vf_2024_07_08.root", "h2_vtx_y_multiplicity" )
    ;
  plot.AddSubPlot( std::vector<double>{ 0.66, 0.0, 0.99, 1.0 } )
    .SetXAxis(Axis().SetTitle("z_{vtx} (mm)").SetLo(-5.0).SetHi(5.0))
    .SetYAxis(Axis().SetTitle("N tracks").SetLo(0).SetHi(500))
    .SetZAxis(Axis().SetTitle("v_{1}").SetLog())
    .AddToPlot( "/home/mikhail/bmn_run8/qa.full_vf_2024_07_08.root", "h2_vtx_z_multiplicity" )
    ;

  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/qa/pictures/vtx_distribution.png" );
};