#include "ris2.h"
#include <Rtypes.h>
#include <TAttMarker.h>
#include <string>

void dd_example(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  std::string file_vf = "~/Flow/BM@N/vf.2024.02.12.root";
  
  auto container= DoubleDifferential{
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
  container.SetStyles( std::vector{
      Style().SetColor( kRed ).SetMarker(kFullCircle),
      Style().SetColor( kBlue+2 ).SetMarker(kFullCircle),
      Style().SetColor( kGreen+2 ).SetMarker(kFullCircle),
      Style().SetColor( kViolet+2 ).SetMarker(kFullCircle),
      Style().SetColor( kCyan+2 ).SetMarker(kFullCircle),
      Style().SetColor( kYellow+2 ).SetMarker(kFullCircle),
  } );

  auto pT_proj = container;
  pT_proj.Rebin( { {"centrality", 1, 10, 30} } )
    .SetSliceAxis({"trProtonY", 6, -0.4, 0.8})
    .SetProjectionAxis({ "trPt", 7, 0.0, 1.4 });
  auto leg_pT = pT_proj.MakeLegend( "y_{cm}"s );
  
  auto y_proj = container;
  y_proj.Rebin( { {"centrality", 1, 10, 30} } )
    .SetSliceAxis({ "trPt", 5, 0.0, 2.0 })
    .SetProjectionAxis({"trProtonY", 10, -0.6, 1.4});
  auto leg_y = y_proj.MakeLegend("p_{T}"s, "GeV/c"s);

  auto plot = Plot( {2000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 0.5, 1.0 } )
      .SetXAxis(Axis().SetTitle("p_{T} (GeV/c)").SetLo(0.0).SetHi(1.5))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.2))
      .AddText( Text().SetText("BM@N RUN8").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddToPlot( *pT_proj )
      .AddLegend( leg_pT );
  plot.AddSubPlot( std::vector<double>{ 0.5, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.6).SetHi(1.4))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.2))
      .AddText( Text().SetText("BM@N RUN8").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddToPlot( *y_proj )
      .AddLegend( leg_y );
  plot.Print( "/home/mikhail/ris2/macro/pictures/dd_example.png" );
}