#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void fit_parameters(){
  using namespace ris2;

  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );

  auto mean_pq_400 = Wrap<Graph>( 
    "mean"s, 
    "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/file_out_400.root"s,
    "g1_mean_2212"s
  );
  mean_pq_400.SetStyle( Style().SetMarker(kFullCircle).SetColor( kBlue ) );

  auto sigma_pq_400 = Wrap<Graph>( 
    "mean"s, 
    "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/file_out_400.root"s,
    "g1_sigma_2212"s
  );
  sigma_pq_400.SetStyle( Style().SetMarker(kFullCircle).SetColor( kBlue ) );

  auto mean_pq_700 = Wrap<Graph>( 
    "mean"s, 
    "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/file_out_700.root"s,
    "g1_mean_2212"s
  );
  mean_pq_700.SetStyle( Style().SetMarker(kFullSquare).SetColor( kRed ) );

  auto sigma_pq_700 = Wrap<Graph>( 
    "mean"s, 
    "/home/mikhail/ris2/macro/bmn_run8/analysis_note/files/file_out_700.root"s,
    "g1_sigma_2212"s
  );
  sigma_pq_700.SetStyle( Style().SetMarker(kFullSquare).SetColor( kRed ) );

  auto leg = new TLegend( 0.2, 0.2, 0.5, 0.5 );
  leg->AddEntry( mean_pq_400.GetResult()->Clone(), "TOF-400", "P" );
  leg->AddEntry( mean_pq_700.GetResult()->Clone(), "TOF-700", "P" );

  auto plot = Plot( {2000, 1100} );

// Upper plots
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 0.5, 1.0 } )
    .SetXAxis(Axis().SetTitle("p/q (GeV/c)").SetLo(0.0).SetHi(6.0))
    .SetYAxis(Axis().SetTitle("#LTm^{2}#GT (GeV^{2}/c^{4})").SetLo(0.6).SetHi(1.0))
    .AddToPlot( mean_pq_400 )
    .AddToPlot( mean_pq_700 )
    .AddLegend(leg)
  ;
  plot.AddSubPlot( std::vector<double>{ 0.5, 0.0, 1.0, 1.0 } )
    .SetXAxis(Axis().SetTitle("p/q (GeV/c)").SetLo(0.0).SetHi(6.0))
    .SetYAxis(Axis().SetTitle("#sigma_{m^{2}} (GeV^{2}/c^{4})").SetLo(-0.0).SetHi(0.5))
    .AddToPlot( sigma_pq_400 )
    .AddToPlot( sigma_pq_700 )
  ;
  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/qa/pictures/fit_parameters.png" );
};