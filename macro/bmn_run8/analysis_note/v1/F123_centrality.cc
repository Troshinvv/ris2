#include "ris2.h"

#include <Rtypes.h>
#include <TAttMarker.h>
#include <string>

void F123_centrality(){
  using namespace ris2;
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );
  
  std::string file_vert = "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/vf.recent.clean.2024.06.18.root";

  auto plot = Plot( {3000, 1100} );

  // **************************
  // F1
  // **************************
  auto R1_f1 = Bunch<Correlation>{};
  R1_f1
    .AddToBunch(
        "F1(F3,T+)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F1_RESCALED(F3_RESCALED,Tpos_RESCALED).y1y1centrality", 
      }
    )
    .AddToBunch(
        "F1(F3,T-)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F1_RESCALED(F3_RESCALED,Tneg_RESCALED).y1y1centrality", 
      }
    )
    .AddToBunch(
        "F1(F2,F3)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality", 
      }
    )
  ;

  R1_f1.GetPalette().SetPalette( std::vector{
      Style().SetColor( kRed ).SetMarker(kFullCircle),
      Style().SetColor( kBlue+2 ).SetMarker(kFullCircle),
      Style().SetColor( kGreen+2 ).SetMarker(kFullCircle),
      Style().SetColor( kRed ).SetMarker(kOpenSquare),
      Style().SetColor( kBlue+2 ).SetMarker(kOpenSquare),
      Style().SetColor( kGreen+2 ).SetMarker(kOpenSquare),
  } );

  auto leg_f1 = R1_f1.MakeLegend( {0.25, 0.85, 0.55, 0.65} );

  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, .33, 1.0 } )
      .SetXAxis(Axis().SetTitle("centrality (%)").SetLo(0.0).SetHi(40.))
      .SetYAxis(Axis().SetTitle("R_{1}").SetLo(0).SetHi(0.55))
      .AddText( Text().SetText("F1").SetPosition({0.25, 0.9}).SetSize(0.05) )
      .AddLegend( leg_f1 )
      .AddToPlot( *R1_f1 );

  // **************************
  // F2
  // **************************
  auto R1_f2 = Bunch<Correlation>{};
  R1_f2
    .AddToBunch(
        "F1{T-}(F1,F3)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F2_RESCALED.Tpos_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality", 
      }
    )
    .AddToBunch(
        "F1{T-}(F1,F3)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F2_RESCALED.Tneg_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality", 
      }
    )
    .AddToBunch(
        "F1(F2,F3)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality", 
      }
    )
  ;
  
  R1_f2.GetPalette().SetPalette( std::vector{
      Style().SetColor( kRed ).SetMarker(kFullCircle),
      Style().SetColor( kBlue+2 ).SetMarker(kFullCircle),
      Style().SetColor( kGreen+2 ).SetMarker(kFullCircle),
      Style().SetColor( kRed ).SetMarker(kOpenSquare),
      Style().SetColor( kBlue+2 ).SetMarker(kOpenSquare),
      Style().SetColor( kGreen+2 ).SetMarker(kOpenSquare),
  } );

  auto leg_f2 = R1_f2.MakeLegend( {0.25, 0.85, 0.55, 0.65} );

  plot.AddSubPlot( std::vector<double>{ 0.33, 0.0, .66, 1.0 } )
      .SetXAxis(Axis().SetTitle("centrality (%)").SetLo(0.0).SetHi(40.))
      .SetYAxis(Axis().SetTitle("R_{1}").SetLo(0).SetHi(0.55))
      .AddText( Text().SetText("F2").SetPosition({0.25, 0.9}).SetSize(0.05) )
      .AddLegend( leg_f2 )
      .AddToPlot( *R1_f2 );

  // **************************
  // F3
  // **************************
  auto R1_f3 = Bunch<Correlation>{};
  R1_f3
    .AddToBunch(
        "F1(F3,T+)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F3_RESCALED(F1_RESCALED,Tpos_RESCALED).y1y1centrality", 
        },
        std::vector<double>{1.}
    )
    .AddToBunch(
        "F1(F3,T-)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F3_RESCALED(F1_RESCALED,Tneg_RESCALED).y1y1centrality", 
      }
    )
    .AddToBunch(
        "F1(F2,F3)"s, 
        file_vert, 
        std::vector<std::string>{ 
          "resolutions/R1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality", 
      }
    )
  ;

  R1_f3.GetPalette().SetPalette( std::vector{
      Style().SetColor( kRed ).SetMarker(kFullCircle),
      Style().SetColor( kBlue+2 ).SetMarker(kFullCircle),
      Style().SetColor( kGreen+2 ).SetMarker(kFullCircle),
      Style().SetColor( kRed ).SetMarker(kOpenSquare),
      Style().SetColor( kBlue+2 ).SetMarker(kOpenSquare),
      Style().SetColor( kGreen+2 ).SetMarker(kOpenSquare),
  } );

  auto leg_f3 = R1_f3.MakeLegend( {0.25, 0.85, 0.55, 0.65} );

  plot.AddSubPlot( std::vector<double>{ 0.66, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("centrality (%)").SetLo(0.0).SetHi(40.))
      .SetYAxis(Axis().SetTitle("R_{1}").SetLo(0).SetHi(0.55))
      .AddText( Text().SetText("F3").SetPosition({0.25, 0.9}).SetSize(0.05) )
      .AddLegend( leg_f2 )
      .AddToPlot( *R1_f3 );

  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/v1/pictures/F123_centrality.png" );
}