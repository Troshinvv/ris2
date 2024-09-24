#include "ris2.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void oven(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  std::string file_vf = "~/Flow/BM@N/vf.2024.02.12.root";
  using namespace ris2;
  
  auto oven_time = std::vector<double>{
    24,
    29,
    34,
    38,
    42,
    47,
    56,
    59,
    2+60,
    5+60,
    7+60,
    9+60,
    12+60,
    15+60,
    17+60,
    20+60 
  };
  auto oven_temp = std::vector<double>{
    26,
    36,
    51,
    63,
    75, 
    88,
    110,
    119,
    127,
    137,
    141,
    149,
    154,
    160,
    164,
    168,
  };
  auto graph_oven = new Graph( oven_time.size(), oven_time.data(), oven_temp.data() );
  graph_oven->Fit("pol1");
  auto graph_data = Wrap<Graph>( graph_oven );
  graph_data.SetStyle( Style().SetColor(kBlack).SetMarker(kFullCircle) );

  auto plot = Plot( {1000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("time (min)").SetLo(20).SetHi(90))
      .SetYAxis(Axis().SetTitle("Temperatume (C)").SetLo(0).SetHi(200))
      .AddToPlot( graph_data )
      ;
  plot.Print( "/home/mikhail/ris2/macro/pictures/oven.png" );
}