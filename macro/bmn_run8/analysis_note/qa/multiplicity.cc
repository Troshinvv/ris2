#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void multiplicity(){
  using namespace ris2;

  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );

  auto with_target = Wrap<Histogram>( 
    "W. target"s, 
    "/home/mikhail/bmn_run8/qa.empty_vf_2024_07_08.root"s,
    "h1_track_multiplicity"s
  );
  with_target.SetStyle( Style().SetMarker(-1).SetColor( kBlue ) );
  with_target.Perform( [](auto p){ p->Scale( 1.0 / p->Integral() ); } );

  auto without_target = Wrap<Histogram>( 
    "W/o target"s, 
    "/home/mikhail/bmn_run8/qa.full_vf_2024_07_08.root"s,
    "h1_track_multiplicity"s
  );
  without_target.SetStyle( Style().SetMarker(-1).SetColor( kRed ) );
  without_target.Perform( [](auto p){ p->Scale( 1.0 / p->Integral() ); } );
  
  without_target.Perform( [](auto p){ 
   auto integral_30pc = p->Integral( 54, 122 ) ;
   auto integral_100pc = p->Integral(1, 250) ;
   std::cout << integral_30pc / integral_100pc << std::endl;
  } );

  auto leg = new TLegend( 0.6, 0.8, 0.9, 0.9 );
  leg->AddEntry( with_target.GetResult()->Clone(), "w. target", "L" );
  leg->AddEntry( without_target.GetResult()->Clone(), "w/o target", "L" );

  auto plot = Plot( {1000, 1100} );

  auto line_30pc = new TLine( 54, .000001, 54, 0.025 );
  auto line_10pc = new TLine( 122, .000001, 122, 0.025 );

// Upper plots
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(Axis().SetTitle("N tracks").SetLo(0.0).SetHi(250.0))
    .SetYAxis(Axis().SetLo(0.000001).SetHi(0.025).SetLog())
    .AddToPlot( with_target )
    .AddToPlot( without_target )
    .AddLegend(leg)
    .AddLine(line_30pc)
    .AddLine(line_10pc)
  ;
  plot.Print( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/qa/pictures/multiplicity.png" );
};