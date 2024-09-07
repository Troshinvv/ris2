#include "ris2.h"
// #include "src/wrap.h"
// #include "src/wrap.h"
#include <Axis.hpp>
#include <Rtypes.h>
#include <TAttMarker.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <string>

void v1_proton_pT_y(){
  using namespace ris2;

  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );
  
  std::vector<double> centrality_edges{ 10, 30 };
  std::vector<double> y_egdes{ 0.5, 1.0 };
  std::vector<double> pT_egdes;

  auto model_y = Wrap<Correlation>{
    "JAM"s,
    "~/Flow/BM@N/jam.2024.02.01.root"s, 
    std::vector<std::string>{
      "tru_proton/v1.psi_rp_PLAIN.x1x1centrality", 
      "tru_proton/v1.psi_rp_PLAIN.y1y1centrality" 
    }
  }.SetStyle( Style().SetColor(kBlack).SetMarker(-1) );
  model_y
    .Rebin( std::vector<Qn::AxisD>{ 
      {"centrality", 1, 10, 30}, 
      {"simPt", 1, 1.0, 2.0} } )
    .Project(std::vector<Qn::AxisD>{{"simProtonY", 8, -0.2, 1.4}});
  
  auto model_pT = Wrap<Correlation>{
    "JAM"s,
    "~/Flow/BM@N/jam.2024.02.01.root"s, 
    std::vector<std::string>{
      "tru_proton/v1.psi_rp_PLAIN.x1x1centrality", 
      "tru_proton/v1.psi_rp_PLAIN.y1y1centrality" 
    }
  }.SetStyle( Style().SetColor(kBlack).SetMarker(-1) );
  model_pT
    .Rebin( std::vector<Qn::AxisD>{ 
      {"centrality", 1, 10, 40},        
      {"simProtonY", 1, 1.0, 1.4}} )
    .Project(std::vector<Qn::AxisD>{{"simPt", 9, 0.2, 2.0},});
  
  std::string file_vf = "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/vf.recent.clean.2024.06.18.root";
  auto v1_data = Wrap<Correlation>{
      "RUN8",
      file_vf,
      std::vector<std::string>{
          // "proton/v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "proton/v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "proton/v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
  v1_data.SetStyle( Style().SetColor(kRed).SetMarker(kFullCircle) );

  auto y_proj = v1_data;
  y_proj
    .Rebin( std::vector<Qn::AxisD>{ 
      {"centrality", 1, 10, 30}, 
      {"trPt", 1, 0.6, 1.4 }, 
    } )
    .Project( std::vector<Qn::AxisD>{{"trProtonY", 8, -0.2, 1.4}} );
  auto fit_function = new TF1("fit1", "pol3", 0.0, 2.0);
  fit_function->FixParameter(2, 0.0);
  y_proj.GetResult()->Print();
  
  auto y_proj_sys = Result<Graph>( dynamic_cast<TGraphErrors*>(y_proj.GetResult()->Clone()) );
  y_proj_sys.Perform( []( auto graph ){
    for( int i=0; i<graph->GetN(); ++i ){
      auto y_err = graph->GetPointY(i) * 0.08;
      graph->SetPointError(i, 0.025, y_err);
    }
  } );
  auto y_graph_sys= Wrap<Graph>();
  y_graph_sys.SetResult( std::move(y_proj_sys) );
  y_graph_sys.SetStyle( Style().SetMarker(kFullCircle).SetColor(kRed) );
  y_graph_sys.SetOption("2"s);

  auto pT_proj = v1_data;
  pT_proj
    .Rebin( std::vector<Qn::AxisD>{ 
      {"centrality", 1, 10, 30}, 
      {"trProtonY", 1, 1.0, 1.4 } 
    } )
    .Project( std::vector<Qn::AxisD>{{"trPt", 9, 0.2, 2.0}} ); 

  auto pT_proj_sys = Result<Graph>( dynamic_cast<TGraphErrors*>(pT_proj.GetResult()->Clone()) );
  pT_proj_sys.Perform( []( auto graph ){
    for( int i=0; i<graph->GetN(); ++i ){
      auto y_err = graph->GetPointY(i) * 0.08;
      graph->SetPointError(i, 0.025, y_err);
    }
  } );
  auto pT_graph_sys= Wrap<Graph>();
  pT_graph_sys.SetResult( std::move(pT_proj_sys) );
  pT_graph_sys.SetStyle( Style().SetMarker(kFullCircle).SetColor(kRed) );
  pT_graph_sys.SetOption("2"s);


  auto star_x = std::vector{
  -0.9502286694874349, 
  -0.7493336863329421, 
  -0.6510663440164145, 
  -0.5516867300666879, 
  -0.45119807313872085, 
  -0.3507029589008367, 
  -0.2502062303354735, 
  -0.15079917281859956,
  -0.05029598694331949, 
  };
  auto star_y = std::vector{
  -0.2376238390967515,
  -0.16518767364109951,
  -0.131034961490218,
  -0.1039736670901586,
  -0.08282210272611482,
  -0.06403391379168397,
  -0.04583656871465633,
  -0.02881961989045176,
  -0.012985650243037009,
  };
  auto graph_star = new Graph( star_x.size(), star_x.data(), star_y.data() );
  auto star_data = Wrap<Graph>( graph_star );
  star_data.ScaleXaxis(-1);
  star_data.Scale(-1);
  star_data.SetStyle( Style().SetColor(kBlue).SetMarker(kFullSquare) );

  auto leg = new TLegend( 0.25, 0.8, 0.55, 0.65 );
  leg->AddEntry(y_proj.GetResult()->Clone(), "BMN Xe+CsI", "P");
  // leg->AddEntry(star_data.GetResult()->Clone(), "STAR-FTX Au+Au", "P");
  leg->AddEntry(model_y.GetResult()->Clone(), "JAM Xe+Cs", "L");
  auto plot = Plot( {2000, 1100} );
  
  auto zero_line_y = new TLine( -0.2, 0.0, 1.4, 0.0 );
  zero_line_y->SetLineStyle(9);

  auto zero_line_pT = new TLine( 0.0, 0.0, 2.0, 0.0 );
  zero_line_pT->SetLineStyle(9);

  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 0.5, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.4))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.05).SetHi(0.6))
      .AddText( Text().SetText("BM@N Preliminary").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddText( Text().SetText("10-30%; 1.0<p_{T}<2.0 (GeV/c)").SetPosition({0.25, 0.85}).SetSize(0.035) )
      .AddToPlot( y_proj )
      .AddToPlot( y_graph_sys )
      .AddToPlot( model_y )
      // .AddToPlot( star_data )
      .AddLegend( leg )
      .AddLine( zero_line_y )
      ;
  
  plot.AddSubPlot( std::vector<double>{ 0.5, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("p_{T} (GeV/c)").SetLo(0.0).SetHi(2.0))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.01).SetHi(0.7))
      .AddText( Text().SetText("BM@N Preliminary").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddText( Text().SetText("10-30%; 1.0<y_{cm}<1.4").SetPosition({0.25, 0.85}).SetSize(0.035) )
      .AddToPlot( pT_proj )
      .AddToPlot( pT_graph_sys )
      .AddLine( zero_line_pT )
      .AddToPlot( model_pT )
      ;
  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/v1/pictures/v1_proton_pT_y.png" );
}