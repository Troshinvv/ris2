#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void multiplicity(){
  using namespace ris2;

  gROOT->Macro( "/home/mikhail/ris2/macro/bmn_run8/analysis_note/style.cc" );

  auto au123 = Wrap<Histogram>( 
    "Au+Au E_{kin}=1.23A GeV"s, 
    "/home/mikhail/ris2/macro/hades/files/AuAu@123GEN9_100k.root"s,
    "event_header_PT3/selected_tof_rpc_hits_PT3"s
  );
  au123.SetStyle( Style().SetMarker(-1).SetColor( kGreen+2 ) );
  // au123.Perform( [](auto p){ p->Scale( 1.0 / p->Integral() ); } );

  auto ag123 = Wrap<Histogram>( 
    "Ag+Ag E_{kin}=1.23A GeV"s, 
    "/home/mikhail/ris2/macro/hades/files/AgAg@123GEN4_100k.root"s,
    "event_header_PT3/selected_tof_rpc_hits_PT3"s
  );
  ag123.SetStyle( Style().SetMarker(-1).SetColor( kRed ) );
  // ag123.Perform( [](auto p){ p->Scale( 1.0 / p->Integral() ); } );

  auto ag158 = Wrap<Histogram>( 
    "Ag+Ag E_{kin}=1.58A GeV"s, 
    "/home/mikhail/ris2/macro/hades/files/AgAg@158GEN4_100k.root"s,
    "event_header_PT3/selected_tof_rpc_hits_PT3"s
  );
  ag158.SetStyle( Style().SetMarker(-1).SetColor( kBlue ) );
  // ag158.Perform( [](auto p){ p->Scale( 1.0 / p->Integral() ); } );

  auto leg = new TLegend( 0.5, 0.8, 0.9, 0.9 );
  leg->AddEntry( au123.GetResult()->Clone(), "Au+Au E_{kin}=1.23A GeV", "L" );
  leg->AddEntry( ag123.GetResult()->Clone(), "Ag+Ag E_{kin}=1.23A GeV", "L" );
  leg->AddEntry( ag158.GetResult()->Clone(), "Ag+Ag E_{kin}=1.58A GeV", "L" );

  auto plot = Plot( {1000, 1100} );

  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
    .SetXAxis(Axis().SetTitle("Hits TOF+RPC").SetLo(0.0).SetHi(300.0))
    .SetYAxis(Axis().SetTitle("counts").SetLo(1.0).SetHi(5000.0).SetLog())
    .AddToPlot( au123 )
    .AddToPlot( ag123 )
    .AddToPlot( ag158 )
    .AddLegend(leg)
  ;
  plot.Print( "/home/mikhail/ris2/macro/hades/qa/pictures/hades_multiplicity_comparison.png" );
};