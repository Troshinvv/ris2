#include "ris2.h"
#include <TAttMarker.h>

void dv1dy_sqrt_snn(){
  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );
  using namespace ris2;

  auto p_prl_5488 = Wrap<Graph>(
    "E895 p_{T}>0 GeV/c: PRL 84 (2000) 5488"s, 
    "/home/mikhail/RISOVALKA/macro/p_PRL_84_2000_5488.root"s,
    "Graph0P"s
  );
  p_prl_5488.SetStyle( Style().SetColor(kViolet+2).SetMarker(-1) );
  // p_prl_5488.GetResult()->SetLineWidth(4);

  // auto p_star_fxt_res = Result<Graph>{};
  double p_star_fxt_x[] = {3.0, 3.2, 3.5, 3.9 };
  double p_star_fxt_y[] = {0.3788184722625617, 0.3045692709846064, 0.2315345144738116,0.16227011327330815};
  double p_star_fxt_x_err[] = {0.0};
  double p_star_fxt_y_err[] = {0.002};
  auto p_star_fxt_res = Result<Graph>( new TGraphErrors( 4, p_star_fxt_x, p_star_fxt_y, nullptr, nullptr) );
  auto p_star_fxt = Wrap<Graph>{};
  p_star_fxt.SetResult( std::move(p_star_fxt_res) );
  p_star_fxt.SetStyle( Style().SetColor(kBlue+2).SetMarker(kFullCircle) );
  p_star_fxt.SetTitle( "STAR-FXT Preliminary 10-40%" );

  double p_fopi_x[] = {2.1, 2.2, 2.24, 2.32, 2.40, 2.51};
  double p_fopi_y[] = {0.4/0.445968, 0.428/0.540316, 0.441/0.614398, 0.428/0.677153, 0.448/0.732502, 0.442/0.804541};
  double p_fopi_x_err[] = {0.0};
  double p_fopi_y_err[] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  auto p_fopi_res = Result<Graph>(new TGraphErrors( 6, p_fopi_x, p_fopi_y, nullptr, p_fopi_y_err ));
  auto p_fopi =  Wrap<Graph>{};
  p_fopi.SetResult( std::move(p_fopi_res) );
  p_fopi.SetStyle( Style().SetColor(kGreen+2).SetMarker(kFullTriangleUp) );

  double p_xe_x[] = { 3.26315};
  // double p_xe_y[] = {2.76506e-01};
  double p_xe_y[] = {2.87418e-01};
  
  double p_xe_x_err[] = {0.025};
  double p_xe_y_err[] = { 2.87418e-01 * 0.08  };

  auto p_xe_res = Result<Graph>( new TGraphErrors( 1, p_xe_x, p_xe_y, p_xe_x_err, p_xe_y_err ) );

  auto p_xe = Wrap<Graph>{};

  p_xe.SetResult( std::move(p_xe_res) );  
  p_xe.SetStyle( Style().SetColor(kRed).SetMarker(kFullSquare) );
  p_xe.SetOption("2");

  auto leg1 = new TLegend(0.45, 0.6, 0.9, 0.9);

  leg1->AddEntry( p_fopi.GetResult(), "FOPI Au+Au midcentral", "P" );
  // leg1->AddEntry( p_prl_5488.GetPoints(), p_prl_5488.GetTitle().c_str(), "L" );
  leg1->AddEntry( p_star_fxt.GetResult(),"STAR-FXT Au+Au 10-40%", "P" );
  leg1->AddEntry( p_xe.GetResult(), "BM@N Xe+CsI 10-30%", "P" );

  auto plot = Plot( {1500, 1100} );
  plot
    .AddSubPlot( std::vector<double>{ 0., 0., 1., 1. } )
    .SetXAxis(Axis().SetTitle("#sqrt{S_{NN}}"))
    .SetYAxis(Axis().SetTitle("dv_{1}/dy|_{y=0}"))
    .AddText( Text().SetText("BM@N Preliminary").SetPosition({0.25, 0.9}).SetSize(0.04) )
    // .AddText( Text().SetText("").SetPosition({0.25, 0.85}).SetSize(0.035) )
    // .AddToPlot( p_prl_5488 )
    .AddToPlot( p_star_fxt )
    .AddToPlot( p_fopi )
    .AddToPlot( p_xe )
    .AddLegend(leg1)
  ;
  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/v1/pictures/dv1dy_sqrt_snn.png" );
}