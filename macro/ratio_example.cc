#include "ris2.h"

#include <Rtypes.h>
#include <TAttMarker.h>

void ratio_example(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  std::string file_vf = "~/Flow/BM@N/vf.2024.02.12.root";
  
  auto container= RatioBuilder<Correlation>{
    std::string{"F1+F2+F3"},
    file_vf, 
    std::vector<std::string>{
          "proton/v1.F1_RESCALED(F3_RESCALED,Tpos_RESCALED).y1y1centrality",
          "proton/v1.F1_RESCALED(F3_RESCALED,Tneg_RESCALED).y1y1centrality",

          "proton/v1.F2_RESCALED.Tpos_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
          "proton/v1.F2_RESCALED.Tneg_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "proton/v1.F3_RESCALED(F1_RESCALED,Tpos_RESCALED).y1y1centrality",
          "proton/v1.F3_RESCALED(F1_RESCALED,Tneg_RESCALED).y1y1centrality",
    }
  };
  container.AddToBunch(
    std::string{"F1"},
    file_vf, 
    std::vector<std::string>{
          "proton/v1.F1_RESCALED(F3_RESCALED,Tpos_RESCALED).y1y1centrality",
          "proton/v1.F1_RESCALED(F3_RESCALED,Tneg_RESCALED).y1y1centrality",
    })
    .AddToBunch(
    std::string{"F2"},
    file_vf, 
    std::vector<std::string>{
          "proton/v1.F2_RESCALED.Tpos_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
          "proton/v1.F2_RESCALED.Tneg_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
    })
    .AddToBunch(
    std::string{"F3"},
    file_vf, 
    std::vector<std::string>{
          "proton/v1.F3_RESCALED(F1_RESCALED,Tpos_RESCALED).y1y1centrality",
          "proton/v1.F3_RESCALED(F1_RESCALED,Tneg_RESCALED).y1y1centrality",
    }
  )
  .SetPalette( 
    Style().SetColor( kBlack ).SetMarker(-1),
    std::vector{
      Style().SetColor( kRed ).SetMarker(kFullCircle),
      Style().SetColor( kBlue+2 ).SetMarker(kFullCircle),
      Style().SetColor( kGreen+2 ).SetMarker(kFullCircle),
  } );
  container( []( Correlation& obj ){
      (*obj).Rebin( { {"centrality", 1, 10, 30}, { "trPt", 1, 0.4, 1.4 } } ).Project({"trProtonY"});
  } );
  auto leg = container.MakeLegend( {0.25, 0.8, 0.55, 0.55} );
  auto plot = Plot({2000, 1100});
  plot.AddRatioPlot( container, {0.0, 0.0, 0.5, 1.0}, {0.5, 0.0, 1.0, 1.0} );
  plot.GetSubPlot( 0 )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.2))
      .AddText( Text().SetText("BM@N RUN8").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddText( Text().SetText("p_{T} > 0.2 GeV/c; 10-30%").SetPosition({0.25, 0.85}).SetSize(0.035) )
      .AddLegend(leg);
  plot.GetSubPlot( 1 )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
      .SetYAxis(Axis().SetTitle("ratio").SetLo(0.85).SetHi(1.15))
      .AddFunction( new TF1("one", "1", -0.2, 1.0) );
  plot.Print( "/home/mikhail/ris2/macro/pictures/ratio_example.png" );
}