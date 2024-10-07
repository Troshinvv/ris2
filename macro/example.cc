#include "ris2.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void example(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  std::string file_vf = "~/Flow/BM@N/vf.2024.02.12.root";
  using namespace ris2;

  auto container = Wrap<Correlation>{
    "proton"s,
    file_vf, 
    std::vector<std::string>{
          "proton/v1.F1_RESCALED(F3_RESCALED,Tpos_RESCALED).y1y1centrality",
          "proton/v1.F1_RESCALED(F3_RESCALED,Tneg_RESCALED).y1y1centrality",

          "proton/v1.F2_RESCALED.Tpos_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
          "proton/v1.F2_RESCALED.Tneg_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "proton/v1.F3_RESCALED(F1_RESCALED,Tpos_RESCALED).y1y1centrality",
          "proton/v1.F3_RESCALED(F1_RESCALED,Tneg_RESCALED).y1y1centrality",
    }, std::vector<double>{1., 1. , 1., 1., 1, 1. }
  };
  container.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlue+2) );
  container.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 30}, { "trPt", 1, 0.2, 1.4 } } )
          .Project(std::vector<Qn::AxisD>{ {"trProtonY", 6, -0.2, 1.0} });

  auto plot = Plot( {1000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.6).SetHi(1.4))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.1))
      .AddToPlot( container )
      .AddSystematics( container )
      ;
  plot.Print( "/home/mikhail/ris2/macro/pictures/example.png" );
}