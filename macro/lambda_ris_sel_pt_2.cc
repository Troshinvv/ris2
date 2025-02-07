#include "ris2.h"
// #include "src/picture.h"
#include <TAttMarker.h>
double massL0 = 1.11568;
double xmin =1.09;
double xmax =1.15;
double background(double *x, double *par) {
    return par[0]*(1. + par[1]*x[0] + par[2]*TMath::Power(x[0],2) + par[3]*TMath::Power(x[0],3) + par[4]*TMath::Power(x[0],4)+ par[5]*TMath::Power(x[0],5));
}

double gaussian(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2)))+par[3]*exp(-0.5*(TMath::Power((x[0]-par[4])/par[5],2)));
}

double fitting_function(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*x[0] + par[5]*TMath::Power(x[0],2) + par[6]*TMath::Power(x[0],3) + par[7]*TMath::Power(x[0],4)+par[8]*TMath::Power(x[0],5));
}
double fitting_function_iter2(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*x[0] + par[5]*TMath::Power(x[0],2) + par[6]*TMath::Power(x[0],3) + par[7]*TMath::Power(x[0],4)+ par[8]*TMath::Power(x[0],5))+par[9]*exp(-0.5*(TMath::Power((x[0]-par[10])/par[11],2)));
}
void lambda_ris_sel_pt_2(){
  gROOT->Macro( "/home/valeriy/ris2/macro/style.cc" );
	 std::string file_vf = "~/bmn_lamb/run8_lambda_correlation_newbins_121024.root";
	 std::string file_vf_0608 = "~/bmn_lamb/run8_lambda_correlation_newbins_0608_121024.root";
	 std::string file_vf_0406 = "~/bmn_lamb/run8_lambda_correlation_newbins_0406_121024.root";
	 std::string file_vf_eff = "~/bmn_lamb/run8_lambda_correlation_newbins_eff_081_161024.root";
         std::string file_vf_0608_eff = "~/bmn_lamb/run8_lambda_correlation_newbins_eff_0608_161024.root";
         std::string file_vf_0406_eff = "~/bmn_lamb/run8_lambda_correlation_newbins_eff_0406_161024.root";

	 std::string file_vf_0412_eff = "~/bmn_lamb/run8_lambda_correlation_newbins_eff_0412_161024.root";
	 std::string file_vf_0612_eff = "~/bmn_lamb/run8_lambda_correlation_newbins_eff_0612_161024_test.root";

	 std::string file_vf_sys_nhits3 = "~/bmn_lamb/run8_lambda_correlation_nhits3_281024.root";
	 std::string file_vf_sys_nhits7 = "~/bmn_lamb/run8_lambda_correlation_nhits7_291024.root";
	 std::string file_vf_sys_softcut = "~/bmn_lamb/run8_lambda_correlation_softcut_061124.root";
         std::string file_vf_sys_hardcut = "~/bmn_lamb/run8_lambda_correlation_hardcut_311024.root";
	 std::string file_vf_sys_runid07550 = "~/bmn_lamb/run8_lambda_correlation_runid07550_061124.root";
	 std::string file_vf_sys_runid7901_all = "~/bmn_lamb/run8_lambda_correlation_runid7901_all_061124.root";
	 std::string file_vf_sys_runid7551_7900 = "~/bmn_lamb/run8_lambda_correlation_runid7551_7900_061124.root";


	 std::string file_vf_sys_NP_nhits7 = "~/bmn_lamb/run8_lambda_correlation_NP_nhits7_240125.root";
	 std::string file_vf_sys_NP_main = "~/bmn_lamb/run8_lambda_correlation_NP_main_270125.root";
	 std::string file_vf_sys_NP_main_fix_recent= "~/bmn_lamb/bmn_run8_lambda_correlation_NP_fix_recentering_290125.root";


	 std::string file_jam_test = "~/bmn_lamb/JAM_lambda_correlation_test_290125.root";


 using namespace ris2; 
Wrap<Correlation> container_eta{
    "lambda_good"s,
    file_vf, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };

Wrap<Correlation> container_NP_eta{
    "lambda_good"s,
    file_vf_sys_NP_main, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_NP_eta_fix_recent{
    "lambda_good"s,
    file_vf_sys_NP_main_fix_recent, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_jam{
    "lambda_good"s,
    file_jam_test, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };

Wrap<Correlation> container_eta_jam_psi{
    "lambda_good"s,
    file_jam_test, std::vector<std::string>{
          "lambda/eta_v1.psi_rp_PLAIN.y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_jam_psi_tru{
    "tru_lambda"s,
    file_jam_test, std::vector<std::string>{
          "lambda/tru_v1.psi_rp_PLAIN.y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_jam_psi_sig{
    "tru_lambda"s,
    file_jam_test, std::vector<std::string>{
          "lambda/v1_sig.psi_rp_PLAIN.y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_jam_psi_bckg{
    "tru_lambda"s,
    file_jam_test, std::vector<std::string>{
          "lambda/v1_bckg.psi_rp_PLAIN.y1y1centrality",
    }
  };

Wrap<Correlation> container_eta_nhits3{
    "lambda_good"s,
    file_vf_sys_nhits3, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_nhits7{
    "lambda_good"s,
    file_vf_sys_nhits7, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_NP_nhits7{
    "lambda_good"s,
    file_vf_sys_NP_nhits7, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_softcut{
    "lambda_good"s,
    file_vf_sys_softcut, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_hardcut{
    "lambda_good"s,
    file_vf_sys_hardcut, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_runid07550{
    "lambda_good"s,
    file_vf_sys_runid07550, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_runid7551_7900{
    "lambda_good"s,
    file_vf_sys_runid7551_7900, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_runid7901_all{
    "lambda_good"s,
    file_vf_sys_runid7901_all, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_0608{
    "lambda_good"s,
    file_vf_0608, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_0406{
    "lambda_good"s,
    file_vf_0406, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };

Wrap<Correlation> container_eta_eff{
    "lambda_good"s,
    file_vf_eff, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_0608_eff{
    "lambda_good"s,
    file_vf_0608_eff, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_0406_eff{
    "lambda_good"s,
    file_vf_0406_eff, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };

Wrap<Correlation> container_eta_0412_eff{
    "lambda_good"s,
    file_vf_0412_eff, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
Wrap<Correlation> container_eta_0612_eff{
    "lambda_good"s,
    file_vf_0612_eff, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
  container_eta.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kBlue+2) );

  container_eta_0608.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlue+2) );
  container_eta_0406.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlue+2) );

  container_eta_eff.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlack) );
  container_eta_0608_eff.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlack) );
  container_eta_0406_eff.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlack) );

  container_eta_0412_eff.SetStyle( Style().SetMarker(kOpenCircle).SetColor(kBlack) );
  container_eta_0612_eff.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlack) );

  container_eta_nhits3.SetStyle( Style().SetMarker(kOpenCircle).SetColor(kRed+2) );
  container_eta_nhits7.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlack) );

  container_eta_softcut.SetStyle( Style().SetMarker(kOpenTriangleUp).SetColor(kGreen+2) );
  container_eta_hardcut.SetStyle( Style().SetMarker(kOpenTriangleDown).SetColor(kGreen-2) );

  container_eta_runid07550.SetStyle( Style().SetMarker(kOpenStar).SetColor(kMagenta+2) );
  container_eta_runid7551_7900.SetStyle( Style().SetMarker(kOpenStar).SetColor(kMagenta-2) );
  container_eta_runid7901_all.SetStyle( Style().SetMarker(kOpenStar).SetColor(kMagenta) );


  container_eta_NP_nhits7.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kBlack) );
  container_NP_eta.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kBlack) );
  container_NP_eta_fix_recent.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kRed) );
  container_eta_jam.SetStyle( Style().SetMarker(kFullCircle).SetColor(kRed) );
  container_eta_jam_psi.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlue) );
  container_eta_jam_psi_tru.SetStyle( Style().SetMarker(-1).SetColor(kBlue) );
  container_eta_jam_psi_sig.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kMagenta) );
  container_eta_jam_psi_bckg.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kGreen) );



double xjam[8]={-0.2509603072983355,-0.08399487836107555,0.08194622279129322,0.24891165172855317,0.4169014084507043,0.582842509603073,0.749807938540333,0.9177976952624841};
double yjam[8]={-0.06562499999999999,-0.015144230769230767,0.020913461538461537,0.06706730769230769,0.10817307692307691,0.15288461538461537,0.2408653846153846,0.34399038461538456};

double ex[10]= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

double xstar[10] ={-0.05059288537549407*(-1),-0.15019762845849802*(-1),-0.24980237154150198*(-1),-0.3509881422924901*(-1),-0.4505928853754941*(-1),-0.5517786561264822*(-1),-0.6513833992094862*(-1),-0.7509881422924901*(-1),-0.850592885375494*(-1),-0.950197628458498*(-1)};
double ystar[10] ={-0.008983451536643027*(-1),-0.033569739952718676*(-1),-0.037825059101654845*(-1),-0.0657210401891253*(-1),-0.08652482269503546*(-1),-0.11725768321513003*(-1),-0.15130023640661938*(-1),-0.19432624113475178*(-1),-0.2401891252955083*(-1),-0.2893617021276596*(-1)};
double eyzero[10]= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

TGraphErrors* v1_bmn_jam = new TGraphErrors(8,xjam,yjam,ex,eyzero);
Wrap<Graph> v1_jam_bmn{"v1_jam_bmn",v1_bmn_jam};
v1_jam_bmn.SetStyle( Style().SetMarker(-1).SetColor(kGreen+2) );

double ex_lamb35[5]= {0.,0.,0.,0.,0.};
double xstar_lamb35[5] ={-0.09594095940959409*(-1),-0.2952029520295203*(-1),-0.49815498154981547*(-1),-0.6937269372693726*(-1),-0.8966789667896679*(-1)};
double ystar_lamb35[5] ={-0.026415094339622643*(-1),-0.05849056603773585*(-1),-0.11132075471698114*(-1),-0.18301886792452832*(-1),-0.2679245283018868*(-1)};
double ey_lamb35[5]= {0.,0.,0.,0.,0.};
TGraphErrors* v1_star_lamb32 = new TGraphErrors(8,xstar,ystar,ex,eyzero);
Wrap<Graph> v1_star32{"v1_star32",v1_star_lamb32};



double ystar081[10] ={-0.01186046511627907*(-1),-0.04534883720930233*(-1),-0.06906976744186047*(-1),-0.10883720930232559*(-1),-0.1458139534883721*(-1),-0.19116279069767442*(-1),-0.24488372093023256*(-1),-0.29651162790697677*(-1),-0.35930232558139535*(-1),-0.4206976744186047*(-1)};
double ystar0608[10] ={-0.010843373493975905*(-1),-0.03493975903614458*(-1),-0.054216867469879526*(-1),-0.08313253012048194*(-1),-0.11445783132530121*(-1),-0.15783132530120483*(-1),-0.1927710843373494*(-1),-0.2421686746987952*(-1),-0.30481927710843376*(-1),-0.35421686746987957*(-1)};
double ystar0406[10] ={-0.00980392156862745*(-1),-0.03431372549019608*(-1),-0.0392156862745098*(-1),-0.06666666666666667*(-1),-0.0892156862745098*(-1),-0.11862745098039215*(-1),-0.15098039215686274*(-1),-0.19509803921568628*(-1),-0.2392156862745098*(-1),-0.2892156862745098*(-1)};

double yBeam_35 = 1.23568;
double yBeam_32 = 1.12741;
double yBeam_326 = 1.15141;
double Tpass_326 = 180.443;
double Tpass_32 = 285.069;
double Tpass_35 = 250.111; 
double ystar35081[5] ={-0.01894736842105263*(-1),-0.05473684210526315*(-1),-0.10105263157894737*(-1),-0.17473684210526316*(-1),-0.27789473684210525*(-1)};
double ystar350608[5] ={-0.014814814814814815*(-1),-0.03851851851851852*(-1),-0.07111111111111112*(-1),-0.14074074074074075*(-1),-0.2251851851851852*(-1)};
double ystar350406[5] ={-0.008823529411764706*(-1),-0.023529411764705882*(-1),-0.055882352941176466*(-1),-0.10735294117647058*(-1),-0.1838235294117647*(-1)};
TGraphErrors* v1_star081_lamb32 = new TGraphErrors(10,xstar,ystar081,ex,eyzero);
Wrap<Graph> v1_star08132{"v1_star08132",v1_star081_lamb32};
TGraphErrors* v1_star0406_lamb32 = new TGraphErrors(10,xstar,ystar0406,ex,eyzero);
Wrap<Graph> v1_star040632{"v1_star040632",v1_star0406_lamb32};
TGraphErrors* v1_star0608_lamb32 = new TGraphErrors(10,xstar,ystar0608,ex,eyzero);
Wrap<Graph> v1_star060832{"v1_star060832",v1_star0608_lamb32};

TGraphErrors* v1_star081_lamb35 = new TGraphErrors(5,xstar_lamb35,ystar35081,ex_lamb35,ey_lamb35);
Wrap<Graph> v1_star08135{"v1_star08135",v1_star081_lamb35};
TGraphErrors* v1_star0406_lamb35 = new TGraphErrors(5,xstar_lamb35,ystar350406,ex_lamb35,ey_lamb35);
Wrap<Graph> v1_star040635{"v1_star040635",v1_star0406_lamb35};
TGraphErrors* v1_star0608_lamb35 = new TGraphErrors(5,xstar_lamb35,ystar350608,ex_lamb35,ey_lamb35);
Wrap<Graph> v1_star060835{"v1_star060835",v1_star0608_lamb35};

v1_star08135.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kGreen+2) );
v1_star060835.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kGreen+2) );
v1_star040635.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kGreen+2) );
v1_star08132.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kRed+2) );
v1_star060832.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kRed+2) );
v1_star040632.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kRed+2) );


TGraphErrors* bmn_081 = container_eta.GetResult();
TGraphErrors* bmn_0608 = container_eta_0608.GetResult();
TGraphErrors* bmn_0406 = container_eta_0406.GetResult();
TGraphErrors* bmn_081_eff = container_eta_eff.GetResult();
TGraphErrors* bmn_0608_eff = container_eta_0608_eff.GetResult();
TGraphErrors* bmn_0406_eff = container_eta_0406_eff.GetResult();
TGraphErrors* bmn_0412_eff = container_eta_0412_eff.GetResult();
TGraphErrors* bmn_0612_eff = container_eta_0612_eff.GetResult();

TGraphErrors* bmn_nhits3 = container_eta_nhits3.GetResult();
TGraphErrors* bmn_nhits7 = container_eta_nhits7.GetResult();
TGraphErrors* bmn_softcut = container_eta_softcut.GetResult();
TGraphErrors* bmn_hardcut = container_eta_hardcut.GetResult();
TGraphErrors* bmn_runid07550 = container_eta_runid07550.GetResult();
TGraphErrors* bmn_runid7551_7900 = container_eta_runid7551_7900.GetResult();
TGraphErrors* bmn_runid7901_all = container_eta_runid7901_all.GetResult();

TGraphErrors* bmn_NP_nhits7 = container_eta_NP_nhits7.GetResult();
TGraphErrors* bmn_NP_eta = container_NP_eta.GetResult();
TGraphErrors* bmn_NP_eta_fix_recent = container_NP_eta_fix_recent.GetResult();

TGraphErrors* bmn_jam_eta = container_eta_jam.GetResult();
TGraphErrors* bmn_jam_eta_psi = container_eta_jam_psi.GetResult();
TGraphErrors* bmn_jam_eta_psi_tru = container_eta_jam_psi_tru.GetResult();

TGraphErrors* star35_081=v1_star08135.GetResult();
TGraphErrors* star35_0608=v1_star060835.GetResult();
TGraphErrors* star35_0406=v1_star040635.GetResult();
TGraphErrors* star32_081=v1_star08132.GetResult();
TGraphErrors* star32_0608=v1_star060832.GetResult();
TGraphErrors* star32_0406=v1_star040632.GetResult();
/*for(int i=0;i<5;i++){
bmn_081->SetPointX(i+1,bmn_081->GetPointX(i+1)/yBeam_326);
bmn_0406->SetPointX(i+1,bmn_0406->GetPointX(i+1)/yBeam_326);
bmn_0608->SetPointX(i+1,bmn_0608->GetPointX(i+1)/yBeam_326);
bmn_081_eff->SetPointX(i+1,bmn_081_eff->GetPointX(i+1)/yBeam_326);
bmn_0406_eff->SetPointX(i+1,bmn_0406_eff->GetPointX(i+1)/yBeam_326);
bmn_0608_eff->SetPointX(i+1,bmn_0608_eff->GetPointX(i+1)/yBeam_326);
bmn_0412_eff->SetPointX(i+1,bmn_0412_eff->GetPointX(i+1)/yBeam_326);
bmn_0612_eff->SetPointX(i+1,bmn_0612_eff->GetPointX(i+1)/yBeam_326);
bmn_nhits3->SetPointX(i+1,bmn_nhits3->GetPointX(i+1)/yBeam_326);
bmn_nhits7->SetPointX(i+1,bmn_nhits7->GetPointX(i+1)/yBeam_326);
bmn_NP_nhits7->SetPointX(i+1,bmn_NP_nhits7->GetPointX(i+1)/yBeam_326);
bmn_NP_eta->SetPointX(i+1,bmn_NP_eta->GetPointX(i+1)/yBeam_326);
bmn_softcut->SetPointX(i+1,bmn_softcut->GetPointX(i+1)/yBeam_326);
bmn_hardcut->SetPointX(i+1,bmn_hardcut->GetPointX(i+1)/yBeam_326);
bmn_runid07550->SetPointX(i+1,bmn_runid07550->GetPointX(i+1)/yBeam_326);
bmn_runid7551_7900->SetPointX(i+1,bmn_runid7551_7900->GetPointX(i+1)/yBeam_326);
bmn_runid7901_all->SetPointX(i+1,bmn_runid7901_all->GetPointX(i+1)/yBeam_326);
star35_081->SetPointX(i+1,star35_081->GetPointX(i+1)/yBeam_35);
star35_0608->SetPointX(i+1,star35_0608->GetPointX(i+1)/yBeam_35);
star35_0406->SetPointX(i+1,star35_0406->GetPointX(i+1)/yBeam_35);

bmn_jam_eta->SetPointX(i+1,bmn_jam_eta->GetPointX(i+1)/yBeam_326);
bmn_jam_eta_psi->SetPointX(i+1,bmn_jam_eta_psi->GetPointX(i+1)/yBeam_326);
bmn_jam_eta_psi_tru->SetPointX(i+1,bmn_jam_eta_psi_tru->GetPointX(i+1)/yBeam_326);
bmn_NP_eta_fix_recent->SetPointX(i+1,bmn_NP_eta_fix_recent->GetPointX(i+1)/yBeam_326);
}
for(int i=0;i<10;i++){
star32_081->SetPointX(i+1,star32_081->GetPointX(i+1)/yBeam_32);
star32_0608->SetPointX(i+1,star32_0608->GetPointX(i+1)/yBeam_32);
star32_0406->SetPointX(i+1,star32_0406->GetPointX(i+1)/yBeam_32);
}*/

char *int_fit_func_dv1dy = new char[100];
TF1* fit_func_dv1dy;
fit_func_dv1dy = new TF1(int_fit_func_dv1dy, "[0]+[1]*x+[2]*x*x*x", -1, 1);
        fit_func_dv1dy->SetNpx(2000);
        std::vector<double> range_vec={-0.2,1};

//container_eta.Fit(fit_func_dv1dy,range_vec);
container_eta_0608.Fit(fit_func_dv1dy,range_vec);
container_eta_0406.Fit(fit_func_dv1dy,range_vec);

//container_eta_eff.Fit(fit_func_dv1dy,range_vec);
//container_eta_0608_eff.Fit(fit_func_dv1dy,range_vec);
//container_eta_0406_eff.Fit(fit_func_dv1dy,range_vec);

container_eta_0412_eff.Fit(fit_func_dv1dy,range_vec);
//container_eta_0612_eff.Fit(fit_func_dv1dy,range_vec);
TGraphErrors* v1_star_lamb35 = new TGraphErrors(5,xstar_lamb35,ystar_lamb35,ex_lamb35,ey_lamb35);
Wrap<Graph> v1_lamb_star35{"v1_lamb_star35",v1_star_lamb35};
v1_lamb_star35.SetStyle( Style().SetMarker(kFullSquare).SetColor(kGreen+1) );
TLegend* leg = new TLegend(0.25,0.7,0.85,0.9);
Text header;
header.SetPosition(std::array{0.25,0.92});
header.SetText("#Lambda, correction on momentum conservation");
header.SetSize(0.027);
Text headerpt;
headerpt.SetPosition(std::array{0.4,0.92});
headerpt.SetText("0.6<p_{T}<1.2 GeV/c");
headerpt.SetSize(0.027);

TFile *inFile_true = new TFile("/home/valeriy/bmn_lamb/JAM_lambda_tru_v1_250924_3.root");
TProfile* tru_v1_prof = (TProfile*)inFile_true->Get("prof1d_v1_y");
double x_true[10];
double y_true[10];
double ex_true[10];
double ey_true[10];
for(int i=1;i<11;i++){
	ex_true[i-1]=0;
	y_true[i-1]=tru_v1_prof->GetBinContent(i);
	ey_true[i-1]=tru_v1_prof->GetBinError(i);
	x_true[i-1]=tru_v1_prof->GetXaxis()->GetBinCenter(i);
}
TGraphErrors* v1_true = new TGraphErrors(10,x_true,y_true,ex_true,ey_true);
Wrap<Graph> v1_true_cont{"v1_true"s,v1_true};
v1_true_cont.SetStyle( Style().SetMarker(-1).SetColor(kBlue+2) );
//  leg->AddEntry(v1_star08132.GetResult(),"STAR, Au+Au, 3.2 GeV","p");
//  leg->AddEntry(v1_star08135.GetResult(),"STAR, Au+Au, 3.5 GeV","p");
//   leg->AddEntry(container_eta_0406_eff.GetResult(),"BM@N, Xe+Cs(I), 3.26 GeV","p");
//  leg->AddEntry(container_eta_0612_eff.GetResult(),"Previous Prod","p");
//  leg->AddEntry(container_eta_nhits3.GetResult(),"nhits>3","p");
//  leg->AddEntry(container_eta_nhits7.GetResult(),"Previous Prod, nhits>7","p");
//  leg->AddEntry(container_NP_eta.GetResult(),"New Prod","p");
//  leg->AddEntry(container_eta_jam.GetResult(),"JAM reco u1Q1","p");
  leg->AddEntry(container_eta_jam_psi.GetResult(),"JAM reco u1#Psi1","p");
  leg->AddEntry(container_eta_jam_psi_tru.GetResult(),"JAM true","l");
  leg->AddEntry(container_eta_jam_psi_sig.GetResult(),"JAM reco signal","p");
  leg->AddEntry(container_eta_jam_psi_bckg.GetResult(),"JAM reco bckg","p");
//  leg->AddEntry(container_NP_eta_fix_recent.GetResult(),"New Prod","p");
//  leg->AddEntry(v1_jam_bmn.GetResult(),"JAM true","l");
//  leg->AddEntry(container_eta_softcut.GetResult(),"soft sel cuts","p");
//  leg->AddEntry(container_eta_hardcut.GetResult(),"hard sel cuts","p");
//  leg->AddEntry(container_eta_runid07550.GetResult(),"RunId<7551","p");
//  leg->AddEntry(container_eta_runid7551_7900.GetResult(),"7551<RunId<7900","p");
//  leg->AddEntry(container_eta_runid7901_all.GetResult(),"7901<RunId","p");

auto plot_eta = Plot( {1000, 1100} );
  plot_eta.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
//	  .AddToPlot( container_eta )
//	  .AddToPlot( container_eta_0406 )
//	  .AddToPlot( container_eta_0608 )
//	  .AddToPlot( container_eta_0406_eff )
//          .AddToPlot( container_eta_0406_eff )
//          .AddToPlot( container_eta_0608_eff )
//	  .AddToPlot( container_eta_0412_eff )
//	  .AddToPlot( container_eta_0612_eff )
//	  .AddToPlot( v1_jam_bmn )
	  .AddToPlot( container_eta_jam_psi )
	  .AddToPlot( container_eta_jam_psi_tru )
	  .AddToPlot( container_eta_jam_psi_sig )
	  .AddToPlot( container_eta_jam_psi_bckg )
//	  .AddToPlot( container_eta_nhits3 )
//	  .AddToPlot( container_eta_nhits7 )
//	  .AddToPlot( container_NP_eta )
//	   .AddToPlot( container_NP_eta_fix_recent )
//         	  .AddToPlot( container_eta_jam )
//	  .AddToPlot( container_eta_softcut )
//	  .AddToPlot( container_eta_hardcut )
//	  .AddToPlot( container_eta_runid07550 )
//	  .AddToPlot( container_eta_runid7551_7900 )
//	  .AddToPlot( container_eta_runid7901_all )
//	  .AddToPlot( v1_true_cont )
//	  .AddToPlot(v1_star08132)
//	  .AddToPlot(v1_star040632)
//	  .AddToPlot(v1_star060832)
//	  .AddToPlot(v1_star040635)
//          .AddToPlot(v1_star08135)
  //        .AddToPlot(v1_star060835)
	  .AddLegend(leg)
	  .AddFunction( new TF1("zero", "0", -0.5, 1.5) )
          .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.2).SetHi(0.5))
 //     .AddSystematics( container_eta )
//      .AddText(header)
      .AddText(headerpt);

  plot_eta.Print( "/home/valeriy/bmn_lamb/hardcut/run8_lambda_v1_eta_comparison_JAM_fix_bckg_peak_310125.png" );
auto container_r= RatioBuilder<Correlation>{
    std::string{"All"},
	    file_vf_0612_eff, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    }
  };
  container_r.AddToBunch(
    std::string{"RunId < 7550"},
    file_vf_sys_runid07550, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    })
    .AddToBunch(
    std::string{"7550 < RunId < 7900"},
    file_vf_sys_runid7551_7900, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    })
   .AddToBunch(
    std::string{"RunId > 7900"},
    file_vf_sys_runid7901_all, std::vector<std::string>{
          "lambda/eta_v1.F1_RESCALED(F2_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F2_RESCALED(F1_RESCALED,F3_RESCALED).y1y1centrality",

          "lambda/eta_v1.F3_RESCALED(F1_RESCALED,F2_RESCALED).y1y1centrality",
    })
  .SetPalette(
    std::vector<Style>{
      Style().SetColor( kBlack ).SetMarker(-1),
      Style().SetColor( kRed ).SetMarker(kFullCircle),
      Style().SetColor( kBlue+2 ).SetMarker(kFullCircle),
      Style().SetColor( kGreen+2 ).SetMarker(kFullCircle),
  } );

auto leg_r = container_r.MakeLegend( {0.25, 0.8, 0.55, 0.55} );
  auto plot_r = Plot({2000, 1100});
  plot_r.AddRatioPlot( container_r, {0.0, 0.0, 0.5, 1.0}, {0.5, 0.0, 1.0, 1.0} );
  plot_r.GetSubPlot( 0 )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
      .SetYAxis(Axis().SetTitle("v_{1}").SetLo(-0.1).SetHi(0.4))
      .AddText( Text().SetText("BM@N RUN8, RunId sys").SetPosition({0.25, 0.9}).SetSize(0.04) )
      .AddText( Text().SetText("0.6 < p_{T} < 1.2 GeV/c; 10-40%").SetPosition({0.25, 0.85}).SetSize(0.035) )
      .AddLegend(leg_r);
  plot_r.GetSubPlot( 1 )
      .SetXAxis(Axis().SetTitle("y_{cm}").SetLo(-0.2).SetHi(1.0))
      .SetYAxis(Axis().SetTitle("ratio").SetLo(0.5).SetHi(2))
      .AddFunction( new TF1("one", "1", -0.2, 1.0) );
  plot_r.Print( "/home/valeriy/bmn_lamb/hardcut/run8_lambda_v1_eta_runid_sys_251224.png" );




double ex_star3[10]= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

double xstar3[10] ={-0.04868913857677903*(-1),-0.149812734082397*(-1),-0.250936329588015*(-1),-0.34831460674157305*(-1),-0.449438202247191*(-1),-0.550561797752809*(-1),-0.6479400749063671*(-1),-0.7490636704119851*(-1),-0.850187265917603*(-1),-0.951310861423221*(-1)};
double ystar3[10] ={-0.024683544303797468*(-1),-0.049367088607594936*(-1),-0.07784810126582278*(-1),-0.110126582278481*(-1),-0.1481012658227848*(-1),-0.19936708860759492*(-1),-0.23924050632911392*(-1),-0.29430379746835444*(-1),-0.3550632911392405*(-1),-0.38164556962025314*(-1)};
double eystar3[10]= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


double ex_lamb39[5]= {0.,0.,0.,0.,0.};
double xstar_lamb39[5] ={-0.09594095940959409*(-1),-0.2952029520295203*(-1),-0.49815498154981547*(-1),-0.6974169741697417*(-1),-0.8966789667896679*(-1)};
double ystar_lamb39[5] ={-0.02088607594936709*(-1),-0.04746835443037974*(-1),-0.07784810126582278*(-1),-0.1310126582278481*(-1),-0.20316455696202532*(-1)};
double ey_lamb39[5]= {0.,0.,0.,0.,0.};
TGraphErrors* v1_star_lamb3 = new TGraphErrors(10,xstar3,ystar3,ex_star3,eystar3);
Wrap<Graph> v1_star3{"v1_star3",v1_star_lamb3};


TGraphErrors* v1_star_lamb39 = new TGraphErrors(5,xstar_lamb39,ystar_lamb39,ex_lamb39,ey_lamb39);
Wrap<Graph> v1_lamb_star39{"v1_lamb_star39",v1_star_lamb39};

double x_bmn_slope[1]={3.26};
double y_bmn_slope[1]={fit_func_dv1dy->GetParameter(1)};
double ex_bmn_slope[1]={0.};
double ey_bmn_slope[1]={fit_func_dv1dy->GetParError(1)};
TGraphErrors* v1_bmn_slope = new TGraphErrors(1,x_bmn_slope,y_bmn_slope,ex_bmn_slope,ey_bmn_slope);
Wrap<Graph> v1_slope_bmn{"v1_slope_bmn"s,v1_bmn_slope};
v1_slope_bmn.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlue+2) );

double y_star_slope_o[4]={0.,0.,0.,0.};
v1_star3.Fit(fit_func_dv1dy,range_vec);
y_star_slope_o[0]=fit_func_dv1dy->GetParameter(1);
v1_star08132.Fit(fit_func_dv1dy,range_vec);
y_star_slope_o[1]=fit_func_dv1dy->GetParameter(1);
v1_star08135.Fit(fit_func_dv1dy,range_vec);
y_star_slope_o[2]=fit_func_dv1dy->GetParameter(1);
v1_lamb_star39.Fit(fit_func_dv1dy,range_vec);
y_star_slope_o[3]=fit_func_dv1dy->GetParameter(1);
double x_star_slope[4]={3,3.2,3.5,3.9};
double y_star_slope[4]={0.310580204778157,0.2416382252559727,0.18703071672354948,0.13447098976109215};
double ex_star_slope[4]={0.,0.,0.,0.};
double ey_star_slope[4]={0.,0.,0.,0.};
TGraphErrors* v1_star_slope_o = new TGraphErrors(4,x_star_slope,y_star_slope_o,ex_star_slope,ey_star_slope);
Wrap<Graph> v1_slope_star_o{"v1_slope_star_o"s,v1_star_slope_o};
v1_slope_star_o.SetStyle( Style().SetMarker(kFullCircle).SetColor(kGreen+2) );
TGraphErrors* v1_star_slope = new TGraphErrors(4,x_star_slope,y_star_slope,ex_star_slope,ey_star_slope);
Wrap<Graph> v1_slope_star{"v1_slope_star"s,v1_star_slope};

TLegend* leg3 = new TLegend(0.4,0.7,0.85,0.9);
//leg3->AddEntry(v1_slope_star.GetResult(),"STAR, Au+Au","p");
  leg3->AddEntry(v1_slope_bmn.GetResult(),"BM@N, Xe+CsI","p");
  leg3->AddEntry(v1_slope_star_o.GetResult(),"STAR, Au+Au","p");

auto plot_dv1dy = Plot( {1000, 1100} );
  plot_dv1dy.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
//	  .AddToPlot( v1_slope_star )
	  .AddToPlot( v1_slope_bmn )
	  .AddToPlot( v1_slope_star_o )
	  .AddLegend(leg3)
          .AddFunction( new TF1("zero", "0", -0.5, 1.5) )
          .SetXAxis(Axis().SetTitle("Collision energy #sqrt{S_{NN}}, GeV").SetLo(2.8).SetHi(4.0))
      .SetYAxis(Axis().SetTitle("dv_{1}/dy|_{y=0}").SetLo(0).SetHi(0.6))
      .AddText(header);

  plot_dv1dy.Print( "/home/valeriy/bmn_lamb/hardcut/run8_lambda_calc2_27kfilesTEST_180924_noeff_fixalg_v1_slope_comparison.png" );


double x_pt_slope[3]={0.5,0.7,0.9};
double y_326_pt_slope[3]={0.,0.,0.};
double ex_pt_slope[3]={0.,0.,0.};
double ey_326_slope[3]={0.,0.,0.};
double y_32_pt_slope[3]={0.,0.,0.};
double ey_32_slope[3]={0.,0.,0.};
double y_35_pt_slope[3]={0.,0.,0.};
double ey_35_slope[3]={0.,0.,0.};
double y_326_eff_pt_slope[3]={0.,0.,0.};
double ey_326_eff_slope[3]={0.,0.,0.};

double x_pt_int_slope[1]={0.6823};
double y_326_int_pt_slope[1]={0.};
double ex_pt_int_slope[1]={0.};
double ey_326_int_slope[1]={0.};
container_eta_0412_eff.Fit(fit_func_dv1dy,range_vec);
y_326_int_pt_slope[0]=fit_func_dv1dy->GetParameter(1);
ey_326_int_slope[0]=fit_func_dv1dy->GetParError(1);

container_eta_0406.Fit(fit_func_dv1dy,range_vec);
y_326_pt_slope[0]=fit_func_dv1dy->GetParameter(1);
ey_326_slope[0]=fit_func_dv1dy->GetParError(1);
container_eta_0608.Fit(fit_func_dv1dy,range_vec);
y_326_pt_slope[1]=fit_func_dv1dy->GetParameter(1);
ey_326_slope[1]=fit_func_dv1dy->GetParError(1);
container_eta.Fit(fit_func_dv1dy,range_vec);
y_326_pt_slope[2]=fit_func_dv1dy->GetParameter(1);
ey_326_slope[2]=fit_func_dv1dy->GetParError(1);

container_eta_0406_eff.Fit(fit_func_dv1dy,range_vec);
y_326_eff_pt_slope[0]=fit_func_dv1dy->GetParameter(1);
ey_326_eff_slope[0]=fit_func_dv1dy->GetParError(1);
container_eta_0608_eff.Fit(fit_func_dv1dy,range_vec);
y_326_eff_pt_slope[1]=fit_func_dv1dy->GetParameter(1);
ey_326_eff_slope[1]=fit_func_dv1dy->GetParError(1);
container_eta_eff.Fit(fit_func_dv1dy,range_vec);
y_326_eff_pt_slope[2]=fit_func_dv1dy->GetParameter(1);
ey_326_eff_slope[2]=fit_func_dv1dy->GetParError(1);

v1_star040632.Fit(fit_func_dv1dy,range_vec);
y_32_pt_slope[0]=fit_func_dv1dy->GetParameter(1);
ey_32_slope[0]=fit_func_dv1dy->GetParError(1);
v1_star060832.Fit(fit_func_dv1dy,range_vec);
y_32_pt_slope[1]=fit_func_dv1dy->GetParameter(1);
ey_32_slope[1]=fit_func_dv1dy->GetParError(1);
v1_star08132.Fit(fit_func_dv1dy,range_vec);
y_32_pt_slope[2]=fit_func_dv1dy->GetParameter(1);
ey_32_slope[2]=fit_func_dv1dy->GetParError(1);

v1_star040635.Fit(fit_func_dv1dy,range_vec);
y_35_pt_slope[0]=fit_func_dv1dy->GetParameter(1);
ey_35_slope[0]=fit_func_dv1dy->GetParError(1);
v1_star060835.Fit(fit_func_dv1dy,range_vec);
y_35_pt_slope[1]=fit_func_dv1dy->GetParameter(1);
ey_35_slope[1]=fit_func_dv1dy->GetParError(1);
v1_star08135.Fit(fit_func_dv1dy,range_vec);
y_35_pt_slope[2]=fit_func_dv1dy->GetParameter(1);
ey_35_slope[2]=fit_func_dv1dy->GetParError(1);

TGraphErrors* v1_bmn_pt_slope = new TGraphErrors(3,x_pt_slope,y_326_pt_slope,ex_pt_slope,ey_326_slope);
Wrap<Graph> v1_slope_pt_bmn{"v1_slope_pt_bmn"s,v1_bmn_pt_slope};
TGraphErrors* v1_star32_pt_slope = new TGraphErrors(3,x_pt_slope,y_32_pt_slope,ex_pt_slope,ey_32_slope);
Wrap<Graph> v1_slope_pt_star32{"v1_slope_pt_star32"s,v1_star32_pt_slope};
TGraphErrors* v1_star35_pt_slope = new TGraphErrors(3,x_pt_slope,y_35_pt_slope,ex_pt_slope,ey_35_slope);
Wrap<Graph> v1_slope_pt_star35{"v1_slope_pt_star35"s,v1_star35_pt_slope};

TGraphErrors* v1_bmn_pt_eff_slope = new TGraphErrors(3,x_pt_slope,y_326_eff_pt_slope,ex_pt_slope,ey_326_eff_slope);
Wrap<Graph> v1_slope_pt_eff_bmn{"v1_slope_pt_eff_bmn"s,v1_bmn_pt_eff_slope};

TGraphErrors* v1_bmn_pt_int_slope = new TGraphErrors(1,x_pt_int_slope,y_326_int_pt_slope,ex_pt_int_slope,ey_326_int_slope);
Wrap<Graph> v1_slope_pt_int_bmn{"v1_slope_pt_int_bmn"s,v1_bmn_pt_int_slope};

v1_slope_pt_bmn.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kBlue+2) );
v1_slope_pt_star32.SetStyle( Style().SetMarker(kFullCircle).SetColor(kBlack) );
v1_slope_pt_star35.SetStyle( Style().SetMarker(kFullTriangleDown).SetColor(kGreen+2) );
v1_slope_pt_eff_bmn.SetStyle( Style().SetMarker(kOpenSquare).SetColor(kRed+2) );
v1_slope_pt_int_bmn.SetStyle( Style().SetMarker(kFullCircle).SetColor(kMagenta) );
TLegend* leg4 = new TLegend(0.2,0.7,0.55,0.9);
  leg4->AddEntry(v1_slope_pt_star32.GetResult(),"STAR, Au+Au,3.2 GeV","p");
  leg4->AddEntry(v1_slope_pt_star35.GetResult(),"STAR, Au+Au,3.5 GeV","p");
  leg4->AddEntry(v1_slope_pt_bmn.GetResult(),"BM@N, Xe+CsI, 3.26 GeV","p");
  leg4->AddEntry(v1_slope_pt_eff_bmn.GetResult(),"BM@N, Xe+CsI, 3.26 GeV,Eff","p");
  leg4->AddEntry(v1_slope_pt_int_bmn.GetResult(),"BM@N, Xe+CsI, 3.26 GeV,0.4<p_{T}<1.2, Eff","p");

auto plot_dv1dy_pt = Plot( {1000, 1100} );
  plot_dv1dy_pt.AddSubPlot( std::vector<double>{ 0.0, 0.0, 1.0, 1.0 } )
          .AddToPlot( v1_slope_pt_star32 )
	  .AddToPlot( v1_slope_pt_star35 )
          .AddToPlot( v1_slope_pt_bmn )
	  .AddToPlot( v1_slope_pt_eff_bmn )
	  .AddToPlot( v1_slope_pt_int_bmn )
          .AddLegend(leg4)
          .AddFunction( new TF1("zero", "0", -0.5, 1.5) )
          .SetXAxis(Axis().SetTitle("p_{T} bins, GeV").SetLo(0.3).SetHi(1.2))
      .SetYAxis(Axis().SetTitle("dv_{1}/dy|_{y=0}").SetLo(-0.3).SetHi(0.6));

  plot_dv1dy_pt.Print( "/home/valeriy/bmn_lamb/hardcut/run8_lambda_calc2_v1_slope_in_diff_pt_bins_comparison_151024.png" );
//  auto points = container.GetResult();
//  TFile FilePoints("/home/valeriy/bmn_lamb/hardcut/pt23_y1_points.root","recreate");
//  points->Write();
//

}
