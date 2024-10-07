#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void v1_fopi_04(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  std::string file_vf = "~/Flow/BM@N/vf.2024.02.12.root";
  using namespace ris2;

  std::vector<double> au_x{0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 };
  std::vector<double> au_y{
    0.018301344610713532,
    0.054215463023205346,
    0.09389909997831275,
    0.13107108002602477,
    0.16698248752982003,
    0.20289253957926703,
    0.23440414226848846,
    0.26339731077857303,
    0.2791829321188463,
    0.28868060073736723,
    0.30572408371286064,
  };

  std::vector<double> ru_x{0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 };
  std::vector<double> ru_y{
    0.0157869767946216,
    0.04666693775753633,
    0.08132183908045981,
    0.11597674040338324,
    0.15125921708956847,
    0.1846549013229235,
    0.21427835610496646,
    0.24138473216222084,
    0.26345830622424643,
    0.2754716981132076,
    0.29188760572543926,
  };

  std::vector<double> ca_x{0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 };
  std::vector<double> ca_y{
    0.012645033615267864,
    0.046039362394274586,
    0.07503253090435918,
    0.10591249186727394,
    0.1349070158317068,
    0.1639001843417914,
    0.1910079158533941,
    0.2111960529169378,
    0.2250948818043809,
    0.2308176100628931,
    0.24534672522229456,

  };
  
  auto graph_au = new Graph( au_x.size(), au_x.data(), au_y.data() );
  auto v1_au = Wrap<Graph>{ "Au+Au", graph_au };
  v1_au.ScaleXaxis( 0.446748 );
  v1_au.SetStyle( Style().SetColor(kRed) );
  v1_au.Fit( new TF1( "au+au", "pol3" ) );

  auto graph_ru = new Graph( ru_x.size(), ru_x.data(), ru_y.data() );
  auto v1_ru = Wrap<Graph>{ "Ru+Ru", graph_ru };
  v1_ru.ScaleXaxis( 0.446748 );
  v1_ru.SetStyle( Style().SetColor(kBlue) );
  v1_ru.Fit( new TF1( "au+au", "pol3" ) );


  auto graph_ca = new Graph( ca_x.size(), ca_x.data(), ca_y.data() );
  auto v1_ca = Wrap<Graph>{ "Ca+Ca", graph_ca };
  v1_ca.ScaleXaxis( 0.446748 );
  v1_ca.SetStyle( Style().SetColor(kGreen+2) );
  v1_ca.Fit( new TF1( "au+au", "pol3" ) );

  auto plot = Plot( {1000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.05).SetHi(1.1))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.5))
      .AddToPlot( v1_au )
      .AddToPlot( v1_ru )
      .AddToPlot( v1_ca )
      ;
  plot.Print( "/home/mikhail/ris2/macro/pictures/v1_fopi_04.png" );
};