#include "ris2.h"

#include <Axis.hpp>
#include <Rtypes.h>
#include <TAttMarker.h>
#include <string>
#include <vector>

void v1_proton_tof(){
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );
  std::string file_100 = "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/vf.recent.2024.06.13.root";
  
  using namespace ris2;

  auto model = Wrap<Correlation>{
    "JAM"s,
    "~/Flow/BM@N/jam.2024.02.01.root"s, std::vector<std::string>{
      "tru_proton/v1.psi_rp_PLAIN.x1x1centrality", 
      "tru_proton/v1.psi_rp_PLAIN.y1y1centrality" 
    }
  }.SetStyle( Style().SetColor(kBlack).SetMarker(-1) );
  model.Rebin( std::vector<Qn::AxisD>{ {"centrality", 1, 10, 30}, 
                    {"simPt", 1, 0.8, 1.6}, } )
        .Project(std::vector<Qn::AxisD>{{"simProtonY", 6, -0.2, 1.0}});
  
  auto vf_whole = Wrap<Correlation>{
    "whole"s,
    file_100, std::vector<std::string>{
      "proton/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",
      // "proton/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
      // "proton/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  }.SetStyle( Style().SetColor(kRed).SetMarker(-2) );
  vf_whole.Rebin( std::vector<Qn::AxisD>{ {"centrality", 1, 10, 30}, {"trPt", 1, 0.8, 1.6}, } )
        .Project(std::vector<Qn::AxisD>{{"trProtonY", 12, -0.2, 1.0}});
  
  auto vf_400 = Wrap<Correlation>{
    "TOF-400"s,
    file_100, std::vector<std::string>{
      "proton_400/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_400/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_400/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  }.SetStyle( Style().SetColor(kBlue).SetMarker(kFullCircle) );
  vf_400.Rebin( std::vector<Qn::AxisD>{ {"centrality", 1, 10, 30}, {"trPt", 1, 0.8, 1.6}, } )
        .Project(std::vector<Qn::AxisD>{{"trProtonY", 9, -0.2, 0.7}});
  

  auto vf_700 = Wrap<Correlation>{
    "TOF-700"s,
    file_100, std::vector<std::string>{
      "proton_700/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_700/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_700/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  }.SetStyle( Style().SetColor(kGreen+2).SetMarker(kFullCircle) );
  vf_700.Rebin( std::vector<Qn::AxisD>{ {"centrality", 1, 10, 30}, {"trPt", 1, 0.8, 1.6}, } )
        .Project(std::vector<Qn::AxisD>{{"trProtonY", 8, 0.2, 1.0}});
  
  auto plot = Plot( {1000, 1100} );

  auto leg = new TLegend( 0.25, 0.8, 0.7, 0.5 );
  // leg->AddEntry( model.GetResult()->Clone(), "JAM", "L" );
  // leg->AddEntry( vf_whole.GetResult()->Clone(), "Whole", "l" );
  leg->AddEntry( vf_400.GetResult()->Clone(), "TOF-400", "P" );
  leg->AddEntry( vf_700.GetResult()->Clone(), "TOF-700", "P" );
  
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.45))
      .AddText( Text().SetText("BM@N RUN8 RESCALED").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddText( Text().SetText("0.8<p_{T}<1.6 (GeV/c); 10-30%").SetPosition({0.25, 0.85}).SetSize(0.035) )
      // .AddToPlot( model )
      // .AddToPlot( vf_whole )
      .AddToPlot( vf_400 )
      .AddToPlot( vf_700 )
      .AddLegend( leg )
      ;

  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/v1/pictures/v1_proton_tof.png" );
}