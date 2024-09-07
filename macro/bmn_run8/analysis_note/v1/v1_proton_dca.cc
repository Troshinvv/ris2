#include "ris2.h"

// #include <Axis.hpp>
#include <Rtypes.h>
#include <TAttMarker.h>
#include <string>
#include <vector>

void v1_proton_dca(){
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );
  auto file_100 = std::string{"/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/vf.recent.dca.2024.06.17.root"};
  
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
  
  auto ratio_builder = RatioBuilder<Correlation>{
    "DCA < 5"s,
    file_100, 
    std::vector<std::string>{
      "proton_dca_5/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_dca_5/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_dca_5/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
  ratio_builder
  .AddToBunch(
    "DCA < 4"s,
    file_100, 
    std::vector<std::string>{
      "proton_dca_4/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_dca_4/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_dca_4/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  )
  .AddToBunch(
    "DCA < 3"s,
    file_100, 
    std::vector<std::string>{
      "proton_dca_3/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_dca_3/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",
      "proton_dca_3/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  )
  ;

  ratio_builder.SetPalette(
    {
      Style().SetColor(kRed).SetMarker(kFullSquare),
      Style().SetColor(kBlue).SetMarker(kFullCircle),
      Style().SetColor(kGreen+2).SetMarker(kFullTriangleDown),
    }
  );

  ratio_builder.Perform([]( auto& obj ){
    obj
    .Rebin( std::vector<Qn::AxisD>{{ "centrality", 1, 10, 30 }, { "trPt", 1, 0.8, 1.6 }} )
    .Project( std::vector<Qn::AxisD>{ { "trProtonY", 12, -0.2, 1.0 } } )
    ;
  });

  auto y_proj_sys = Result<Graph>( dynamic_cast<TGraphErrors*>( ratio_builder.GetResultWraps().at(0).GetResult()->Clone() ) );
  y_proj_sys.Perform( []( auto graph ){
    for( int i=0; i<graph->GetN(); ++i ){
      auto y_err = graph->GetPointY(i) * 0.02;
      graph->SetPointError(i, 0.025, y_err);
    }
  } );
  auto y_graph_sys= Wrap<Graph>();
  y_graph_sys.SetResult( std::move(y_proj_sys) );
  y_graph_sys.SetStyle( Style().SetMarker(-1).SetColor(kRed) );
  y_graph_sys.SetOption("3"s);

  
  auto plot = Plot( {1000, 1100} );

  auto leg =  ratio_builder.MakeLegend( {0.25, 0.8, 0.7, 0.5} );
  
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
    .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.45))
    .AddText( Text().SetText("BM@N RUN8").SetPosition({0.25, 0.9}).SetSize(0.04) )
    .AddText( Text().SetText("0.8<p_{T}<1.6 (GeV/c); 10-30%").SetPosition({0.25, 0.85}).SetSize(0.035) )
    .AddToPlot(y_graph_sys)
    .AddToPlot(ratio_builder.GetResultWraps())
    .AddLegend( leg )
  ;
  // plot
  //   .GetSubPlot(1)
  //   .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
  //   .SetYAxis(Axis().SetTitle("ratio").SetLo(0.9).SetHi(1.1))
  // ;

  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/v1/pictures/v1_proton_dca.png" );
}