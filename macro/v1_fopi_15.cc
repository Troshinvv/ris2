#include "ris2_noqntools.h"
// #include "src/picture.h"
#include <TAttMarker.h>

void v1_fopi(){
  gROOT->Macro( "/home/mikhail/ris2/macro/style.cc" );
  std::string file_vf = "~/Flow/BM@N/vf.2024.02.12.root";
  using namespace ris2;

  std::vector<double> au_x{0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 };
  std::vector<double> au_y{
    0.018648018648018683,
    0.05944055944055948, 
    0.10023310023310028,
    0.14452214452214457, 
    0.19114219114219116, 
    0.23776223776223782,
    0.28438228438228447, 
    0.3275058275058276, 
    0.3682983682983684,
    0.41025641025641035,
     0.4440559440559441
  };

  std::vector<double> ru_x{0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 };
  std::vector<double> ru_y{
    0.011655011655011704,
    0.044289044289044344, 
    0.0815850815850816,
    0.12237762237762245,
    0.16433566433566438,
    0.21095571095571097,
    0.2552447552447553,
    0.29370629370629375,
    0.33216783216783224,
    0.36480186480186483,
    0.39393939393939403,
  };

  std::vector<double> ca_x{0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 };
  std::vector<double> ca_y{
    0.010489510489510523,
    0.03962703962703962,
    0.07226107226107231,
    0.10722610722610726,
    0.14219114219114226,
    0.1783216783216784,
    0.2121212121212122,
    0.24358974358974367,
    0.26689976689976697,
    0.2890442890442891,
    0.3076923076923077,
  };
  
  auto graph_au = new Graph( au_x.size(), au_x.data(), au_y.data() );
  auto v1_au = Wrap<Graph>{ "Au+Au", graph_au };
  v1_au.ScaleXaxis( 0.80 );
  v1_au.SetStyle( Style().SetColor(kRed) );
  v1_au.Fit( new TF1( "au+au", "pol1" ) );

  auto graph_ru = new Graph( ru_x.size(), ru_x.data(), ru_y.data() );
  auto v1_ru = Wrap<Graph>{ "Ru+Ru", graph_ru };
  v1_ru.ScaleXaxis( 0.80 );
  v1_ru.SetStyle( Style().SetColor(kBlue) );
  v1_ru.Fit( new TF1( "au+au", "pol1" ) );


  auto graph_ca = new Graph( ca_x.size(), ca_x.data(), ca_y.data() );
  auto v1_ca = Wrap<Graph>{ "Ca+Ca", graph_ca };
  v1_ca.ScaleXaxis( 0.80 );
  v1_ca.SetStyle( Style().SetColor(kGreen+2) );
  v1_ca.Fit( new TF1( "au+au", "pol1" ) );

  auto plot = Plot( {1000, 1100} );
  plot.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.05).SetHi(1.1))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.5))
      .AddToPlot( v1_au )
      .AddToPlot( v1_ru )
      .AddToPlot( v1_ca )
      ;
  plot.Print( "/home/mikhail/ris2/macro/pictures/v1_fopi.png" );
};