
#include "ConfigParser.h"
#include "ParserUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THStack.h"
#include "TPad.h"
#include "TEfficiency.h"
#include "TTree.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

bool passCutBasedJetId(const float& pt,const float& eta, const float& betastarclassic, const float& dR2Mean, const int& nvtx);
std::vector<int> chooseMaxptJets(std::map<int,TLorentzVector>* JetP4);
std::vector<int> chooseMaxptCouple(std::map<int,TLorentzVector>* JetP4);
std::vector<int> chooseMaxptOverMCouple(std::map<int,TLorentzVector>* JetP4);
std::vector<int> chooseMinDeltaMCouple(std::map<int,TLorentzVector>* JetP4,const float mgg);

int main(int argc, char** argv)
{

  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc != 5)
  {
    std::cerr << ">>>>> Analyses::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputFile = gConfigParser -> readStringOption("Input::inputFile");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");

  float dR = atof(argv[2]); 
  int isBtagged = atoi(argv[3]);
  int SelectJet = atoi(argv[4]); 

  std::cout << "inputFile = " << inputFile << std::endl; 
  std::cout << "inputTree = " << inputTree << std::endl; 
  std::cout << "dR        = " << dR << std::endl; 
  std::cout << "isBtagged = " << isBtagged << std::endl; 
  std::cout << "SelectJet = " << SelectJet << std::endl; 
  
  TTree* ntu;
  TFile* file = TFile::Open(inputFile.c_str());
  ntu = (TTree*)file->Get(inputTree.c_str()); 
 
  if(ntu->GetEntries() == 0 )
  {
      std::cout << "Error: input file is empty" << std::endl; 
      return -1;
  }

  // Declaration of leaf types
  int           itype;
  int           run;
  int           lumis;
  int           event;
  float         weight;
  float         evweight;
  float         pu_weight;
  float         pu_n;
  float         nvtx;
  float         rho;
  int           category;
  float         ph1_e;
  float         ph2_e;
  float         ph1_pt;
  float         ph2_pt;
  float         ph1_phi;
  float         ph2_phi;
  float         ph1_eta;
  float         ph2_eta;
  float         j1_e,j2_e,j3_e,j4_e,j5_e,j6_e,j7_e,j8_e,j9_e,j10_e,j11_e,j12_e,j13_e,j14_e,j15_e;
  float         j1_pt,j2_pt,j3_pt,j4_pt,j5_pt,j6_pt,j7_pt,j8_pt,j9_pt,j10_pt,j11_pt,j12_pt,j13_pt,j14_pt,j15_pt;
  float         j1_eta,j2_eta,j3_eta,j4_eta,j5_eta,j6_eta,j7_eta,j8_eta,j9_eta,j10_eta,j11_eta,j12_eta,j13_eta,j14_eta,j15_eta;
  float         j1_phi,j2_phi,j3_phi,j4_phi,j5_phi,j6_phi,j7_phi,j8_phi,j9_phi,j10_phi,j11_phi,j12_phi,j13_phi,j14_phi,j15_phi;
  float         j1_csvBtag,j2_csvBtag,j3_csvBtag,j4_csvBtag,j5_csvBtag,j6_csvBtag,j7_csvBtag,j8_csvBtag,j9_csvBtag,j10_csvBtag,j11_csvBtag, j12_csvBtag,j13_csvBtag,j14_csvBtag,j15_csvBtag;
  float         j1_betaStarClassic,j2_betaStarClassic,j3_betaStarClassic,j4_betaStarClassic,j5_betaStarClassic,j6_betaStarClassic, j7_betaStarClassic,j8_betaStarClassic,j9_betaStarClassic,j10_betaStarClassic,j11_betaStarClassic,j12_betaStarClassic,j13_betaStarClassic, j14_betaStarClassic,j15_betaStarClassic;
  float         j1_dR2Mean,j2_dR2Mean,j3_dR2Mean,j4_dR2Mean,j5_dR2Mean,j6_dR2Mean,j7_dR2Mean,j8_dR2Mean,j9_dR2Mean,j10_dR2Mean,j11_dR2Mean, j12_dR2Mean,j13_dR2Mean,j14_dR2Mean,j15_dR2Mean;
  float         gr_b1_p4_pt;
  float         gr_b1_p4_eta;
  float         gr_b1_p4_phi;
  float         gr_b1_p4_mass;
  float         gr_b2_p4_pt;
  float         gr_b2_p4_eta;
  float         gr_b2_p4_phi;
  float         gr_b2_p4_mass;
  float         gr_hbb_p4_pt;
  float         gr_hbb_p4_eta;
  float         gr_hbb_p4_phi;
  float         gr_hbb_p4_mass;

  // List of branches
  TBranch        *b_itype;   //!
  TBranch        *b_run;   //!
  TBranch        *b_lumis;   //!
  TBranch        *b_event;   //!
  TBranch        *b_weight;   //!
  TBranch        *b_evweight;   //!
  TBranch        *b_pu_weight;   //!
  TBranch        *b_pu_n;   //!
  TBranch        *b_nvtx;   //!
  TBranch        *b_rho;   //!
  TBranch        *b_category;   //!
  TBranch        *b_ph1_e;   //!
  TBranch        *b_ph2_e;   //!
  TBranch        *b_ph1_pt;   //!
  TBranch        *b_ph2_pt;   //!
  TBranch        *b_ph1_phi;   //!
  TBranch        *b_ph2_phi;   //!
  TBranch        *b_ph1_eta;   //!
  TBranch        *b_ph2_eta;   //!
  TBranch        *b_j1_e,*b_j2_e,*b_j3_e,*b_j4_e,*b_j5_e,*b_j6_e,*b_j7_e,*b_j8_e,*b_j9_e,*b_j10_e,*b_j11_e,*b_j12_e,*b_j13_e,*b_j14_e, *b_j15_e;   //!
  TBranch        *b_j1_pt,*b_j2_pt,*b_j3_pt,*b_j4_pt,*b_j5_pt,*b_j6_pt,*b_j7_pt,*b_j8_pt,*b_j9_pt,*b_j10_pt,*b_j11_pt,*b_j12_pt,*b_j13_pt, *b_j14_pt,*b_j15_pt;   //!
  TBranch        *b_j1_eta,*b_j2_eta,*b_j3_eta,*b_j4_eta,*b_j5_eta,*b_j6_eta,*b_j7_eta,*b_j8_eta,*b_j9_eta,*b_j10_eta,*b_j11_eta,*b_j12_eta, *b_j13_eta,*b_j14_eta,*b_j15_eta;   //!
  TBranch        *b_j1_phi,*b_j2_phi,*b_j3_phi,*b_j4_phi,*b_j5_phi,*b_j6_phi,*b_j7_phi,*b_j8_phi,*b_j9_phi,*b_j10_phi,*b_j11_phi,*b_j12_phi, *b_j13_phi,*b_j14_phi,*b_j15_phi;   //!
  TBranch        *b_j1_csvBtag,*b_j2_csvBtag,*b_j3_csvBtag,*b_j4_csvBtag,*b_j5_csvBtag,*b_j6_csvBtag,*b_j7_csvBtag,*b_j8_csvBtag,*b_j9_csvBtag, *b_j10_csvBtag,*b_j11_csvBtag,*b_j12_csvBtag,*b_j13_csvBtag,*b_j14_csvBtag, *b_j15_csvBtag;   //!
  TBranch        *b_j1_betaStarClassic,*b_j2_betaStarClassic,*b_j3_betaStarClassic,*b_j4_betaStarClassic,*b_j5_betaStarClassic, *b_j6_betaStarClassic,*b_j7_betaStarClassic,*b_j8_betaStarClassic,*b_j9_betaStarClassic,*b_j10_betaStarClassic,*b_j11_betaStarClassic, *b_j12_betaStarClassic, *b_j13_betaStarClassic,*b_j14_betaStarClassic, *b_j15_betaStarClassic;   //!
  TBranch        *b_j1_dR2Mean,*b_j2_dR2Mean,*b_j3_dR2Mean,*b_j4_dR2Mean,*b_j5_dR2Mean,*b_j6_dR2Mean,*b_j7_dR2Mean,*b_j8_dR2Mean, *b_j9_dR2Mean,*b_j10_dR2Mean,*b_j11_dR2Mean,*b_j12_dR2Mean,*b_j13_dR2Mean,*b_j14_dR2Mean, *b_j15_dR2Mean;   //!
  TBranch        *b_gr_b1_p4_pt;   //!
  TBranch        *b_gr_b1_p4_eta;   //!
  TBranch        *b_gr_b1_p4_phi;   //!
  TBranch        *b_gr_b1_p4_mass;   //!
  TBranch        *b_gr_b2_p4_pt;   //!
  TBranch        *b_gr_b2_p4_eta;   //!
  TBranch        *b_gr_b2_p4_phi;   //!
  TBranch        *b_gr_b2_p4_mass;   //!
  TBranch        *b_gr_hbb_p4_pt;   //!
  TBranch        *b_gr_hbb_p4_eta;   //!
  TBranch        *b_gr_hbb_p4_phi;   //!
  TBranch        *b_gr_hbb_p4_mass;   //!

  // Set branch addresses and branch pointers
  ntu->SetBranchAddress("itype", &itype, &b_itype);
  ntu->SetBranchAddress("run", &run, &b_run);
  ntu->SetBranchAddress("lumis", &lumis, &b_lumis);
  ntu->SetBranchAddress("event", &event, &b_event);
  ntu->SetBranchAddress("weight", &weight, &b_weight);
  ntu->SetBranchAddress("evweight", &evweight, &b_evweight);
  ntu->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  ntu->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
  ntu->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  ntu->SetBranchAddress("rho", &rho, &b_rho);
  ntu->SetBranchAddress("category", &category, &b_category);
  ntu->SetBranchAddress("ph1_e", &ph1_e, &b_ph1_e);
  ntu->SetBranchAddress("ph2_e", &ph2_e, &b_ph2_e);
  ntu->SetBranchAddress("ph1_pt", &ph1_pt, &b_ph1_pt);
  ntu->SetBranchAddress("ph2_pt", &ph2_pt, &b_ph2_pt);
  ntu->SetBranchAddress("ph1_phi", &ph1_phi, &b_ph1_phi);
  ntu->SetBranchAddress("ph2_phi", &ph2_phi, &b_ph2_phi);
  ntu->SetBranchAddress("ph1_eta", &ph1_eta, &b_ph1_eta);
  ntu->SetBranchAddress("ph2_eta", &ph2_eta, &b_ph2_eta);
  ntu->SetBranchAddress("j1_e", &j1_e, &b_j1_e);
  ntu->SetBranchAddress("j1_pt", &j1_pt, &b_j1_pt);
  ntu->SetBranchAddress("j1_phi", &j1_phi, &b_j1_phi);
  ntu->SetBranchAddress("j1_eta", &j1_eta, &b_j1_eta);
  ntu->SetBranchAddress("j1_csvBtag", &j1_csvBtag, &b_j1_csvBtag);
  ntu->SetBranchAddress("j1_betaStarClassic", &j1_betaStarClassic, &b_j1_betaStarClassic);
  ntu->SetBranchAddress("j1_dR2Mean", &j1_dR2Mean, &b_j1_dR2Mean);
  ntu->SetBranchAddress("j2_e", &j2_e, &b_j2_e);
  ntu->SetBranchAddress("j2_pt", &j2_pt, &b_j2_pt);
  ntu->SetBranchAddress("j2_phi", &j2_phi, &b_j2_phi);
  ntu->SetBranchAddress("j2_eta", &j2_eta, &b_j2_eta);
  ntu->SetBranchAddress("j2_csvBtag", &j2_csvBtag, &b_j2_csvBtag);
  ntu->SetBranchAddress("j2_betaStarClassic", &j2_betaStarClassic, &b_j2_betaStarClassic);
  ntu->SetBranchAddress("j2_dR2Mean", &j2_dR2Mean, &b_j2_dR2Mean);
  ntu->SetBranchAddress("j3_e", &j3_e, &b_j3_e);
  ntu->SetBranchAddress("j3_pt", &j3_pt, &b_j3_pt);
  ntu->SetBranchAddress("j3_phi", &j3_phi, &b_j3_phi);
  ntu->SetBranchAddress("j3_eta", &j3_eta, &b_j3_eta);
  ntu->SetBranchAddress("j3_csvBtag", &j3_csvBtag, &b_j3_csvBtag);
  ntu->SetBranchAddress("j3_betaStarClassic", &j3_betaStarClassic, &b_j3_betaStarClassic);
  ntu->SetBranchAddress("j3_dR2Mean", &j3_dR2Mean, &b_j3_dR2Mean);
  ntu->SetBranchAddress("j4_e", &j4_e, &b_j4_e);
  ntu->SetBranchAddress("j4_pt", &j4_pt, &b_j4_pt);
  ntu->SetBranchAddress("j4_phi", &j4_phi, &b_j4_phi);
  ntu->SetBranchAddress("j4_eta", &j4_eta, &b_j4_eta);
  ntu->SetBranchAddress("j4_csvBtag", &j4_csvBtag, &b_j4_csvBtag);
  ntu->SetBranchAddress("j4_betaStarClassic", &j4_betaStarClassic, &b_j4_betaStarClassic);
  ntu->SetBranchAddress("j4_dR2Mean", &j4_dR2Mean, &b_j4_dR2Mean);
  ntu->SetBranchAddress("j5_e", &j5_e, &b_j5_e);
  ntu->SetBranchAddress("j5_pt", &j5_pt, &b_j5_pt);
  ntu->SetBranchAddress("j5_phi", &j5_phi, &b_j5_phi);
  ntu->SetBranchAddress("j5_eta", &j5_eta, &b_j5_eta);
  ntu->SetBranchAddress("j5_csvBtag", &j5_csvBtag, &b_j5_csvBtag);
  ntu->SetBranchAddress("j5_betaStarClassic", &j5_betaStarClassic, &b_j5_betaStarClassic);
  ntu->SetBranchAddress("j5_dR2Mean", &j5_dR2Mean, &b_j5_dR2Mean);
  ntu->SetBranchAddress("j6_e", &j6_e, &b_j6_e);
  ntu->SetBranchAddress("j6_pt", &j6_pt, &b_j6_pt);
  ntu->SetBranchAddress("j6_phi", &j6_phi, &b_j6_phi);
  ntu->SetBranchAddress("j6_eta", &j6_eta, &b_j6_eta);
  ntu->SetBranchAddress("j6_csvBtag", &j6_csvBtag, &b_j6_csvBtag);
  ntu->SetBranchAddress("j6_betaStarClassic", &j6_betaStarClassic, &b_j6_betaStarClassic);
  ntu->SetBranchAddress("j6_dR2Mean", &j6_dR2Mean, &b_j6_dR2Mean);
  ntu->SetBranchAddress("j7_e", &j7_e, &b_j7_e);
  ntu->SetBranchAddress("j7_pt", &j7_pt, &b_j7_pt);
  ntu->SetBranchAddress("j7_phi", &j7_phi, &b_j7_phi);
  ntu->SetBranchAddress("j7_eta", &j7_eta, &b_j7_eta);
  ntu->SetBranchAddress("j7_csvBtag", &j7_csvBtag, &b_j7_csvBtag);
  ntu->SetBranchAddress("j7_betaStarClassic", &j7_betaStarClassic, &b_j7_betaStarClassic);
  ntu->SetBranchAddress("j7_dR2Mean", &j7_dR2Mean, &b_j7_dR2Mean);
  ntu->SetBranchAddress("j8_e", &j8_e, &b_j8_e);
  ntu->SetBranchAddress("j8_pt", &j8_pt, &b_j8_pt);
  ntu->SetBranchAddress("j8_phi", &j8_phi, &b_j8_phi);
  ntu->SetBranchAddress("j8_eta", &j8_eta, &b_j8_eta);
  ntu->SetBranchAddress("j8_csvBtag", &j8_csvBtag, &b_j8_csvBtag);
  ntu->SetBranchAddress("j8_betaStarClassic", &j8_betaStarClassic, &b_j8_betaStarClassic);
  ntu->SetBranchAddress("j8_dR2Mean", &j8_dR2Mean, &b_j8_dR2Mean);
  ntu->SetBranchAddress("j9_e", &j9_e, &b_j9_e);
  ntu->SetBranchAddress("j9_pt", &j9_pt, &b_j9_pt);
  ntu->SetBranchAddress("j9_phi", &j9_phi, &b_j9_phi);
  ntu->SetBranchAddress("j9_eta", &j9_eta, &b_j9_eta);
  ntu->SetBranchAddress("j9_csvBtag", &j9_csvBtag, &b_j9_csvBtag);
  ntu->SetBranchAddress("j9_betaStarClassic", &j9_betaStarClassic, &b_j9_betaStarClassic);
  ntu->SetBranchAddress("j9_dR2Mean", &j9_dR2Mean, &b_j9_dR2Mean);
  ntu->SetBranchAddress("j10_e", &j10_e, &b_j10_e);
  ntu->SetBranchAddress("j10_pt", &j10_pt, &b_j10_pt);
  ntu->SetBranchAddress("j10_phi", &j10_phi, &b_j10_phi);
  ntu->SetBranchAddress("j10_eta", &j10_eta, &b_j10_eta);
  ntu->SetBranchAddress("j10_csvBtag", &j10_csvBtag, &b_j10_csvBtag);
  ntu->SetBranchAddress("j10_betaStarClassic", &j10_betaStarClassic, &b_j10_betaStarClassic);
  ntu->SetBranchAddress("j10_dR2Mean", &j10_dR2Mean, &b_j10_dR2Mean);
  ntu->SetBranchAddress("j11_e", &j11_e, &b_j11_e);
  ntu->SetBranchAddress("j11_pt", &j11_pt, &b_j11_pt);
  ntu->SetBranchAddress("j11_phi", &j11_phi, &b_j11_phi);
  ntu->SetBranchAddress("j11_eta", &j11_eta, &b_j11_eta);
  ntu->SetBranchAddress("j11_csvBtag", &j11_csvBtag, &b_j11_csvBtag);  
  ntu->SetBranchAddress("j11_betaStarClassic", &j11_betaStarClassic, &b_j11_betaStarClassic);
  ntu->SetBranchAddress("j11_dR2Mean", &j11_dR2Mean, &b_j11_dR2Mean);
  ntu->SetBranchAddress("j12_e", &j12_e, &b_j12_e);
  ntu->SetBranchAddress("j12_pt", &j12_pt, &b_j12_pt);
  ntu->SetBranchAddress("j12_phi", &j12_phi, &b_j12_phi);
  ntu->SetBranchAddress("j12_eta", &j12_eta, &b_j12_eta);
  ntu->SetBranchAddress("j12_csvBtag", &j12_csvBtag, &b_j12_csvBtag);
  ntu->SetBranchAddress("j12_betaStarClassic", &j12_betaStarClassic, &b_j12_betaStarClassic);
  ntu->SetBranchAddress("j12_dR2Mean", &j12_dR2Mean, &b_j12_dR2Mean);
  ntu->SetBranchAddress("j13_e", &j13_e, &b_j13_e);
  ntu->SetBranchAddress("j13_pt", &j13_pt, &b_j13_pt);
  ntu->SetBranchAddress("j13_phi", &j13_phi, &b_j13_phi);
  ntu->SetBranchAddress("j13_eta", &j13_eta, &b_j13_eta);
  ntu->SetBranchAddress("j13_csvBtag", &j13_csvBtag, &b_j13_csvBtag);
  ntu->SetBranchAddress("j13_betaStarClassic", &j13_betaStarClassic, &b_j13_betaStarClassic);
  ntu->SetBranchAddress("j13_dR2Mean", &j13_dR2Mean, &b_j13_dR2Mean);
  ntu->SetBranchAddress("j14_e", &j14_e, &b_j14_e);
  ntu->SetBranchAddress("j14_pt", &j14_pt, &b_j14_pt);
  ntu->SetBranchAddress("j14_phi", &j14_phi, &b_j14_phi);
  ntu->SetBranchAddress("j14_eta", &j14_eta, &b_j14_eta);
  ntu->SetBranchAddress("j14_csvBtag", &j14_csvBtag, &b_j14_csvBtag);
  ntu->SetBranchAddress("j14_betaStarClassic", &j14_betaStarClassic, &b_j14_betaStarClassic);
  ntu->SetBranchAddress("j14_dR2Mean", &j14_dR2Mean, &b_j14_dR2Mean);
  ntu->SetBranchAddress("j15_e", &j15_e, &b_j15_e);
  ntu->SetBranchAddress("j15_pt", &j15_pt, &b_j15_pt);
  ntu->SetBranchAddress("j15_phi", &j15_phi, &b_j15_phi);
  ntu->SetBranchAddress("j15_eta", &j15_eta, &b_j15_eta);
  ntu->SetBranchAddress("j15_csvBtag", &j15_csvBtag, &b_j15_csvBtag);
  ntu->SetBranchAddress("j15_betaStarClassic", &j15_betaStarClassic, &b_j15_betaStarClassic);
  ntu->SetBranchAddress("j15_dR2Mean", &j15_dR2Mean, &b_j15_dR2Mean);
  ntu->SetBranchAddress("gr_b1_p4_pt", &gr_b1_p4_pt, &b_gr_b1_p4_pt);
  ntu->SetBranchAddress("gr_b1_p4_eta", &gr_b1_p4_eta, &b_gr_b1_p4_eta);
  ntu->SetBranchAddress("gr_b1_p4_phi", &gr_b1_p4_phi, &b_gr_b1_p4_phi);
  ntu->SetBranchAddress("gr_b1_p4_mass", &gr_b1_p4_mass, &b_gr_b1_p4_mass);
  ntu->SetBranchAddress("gr_b2_p4_pt", &gr_b2_p4_pt, &b_gr_b2_p4_pt);
  ntu->SetBranchAddress("gr_b2_p4_eta", &gr_b2_p4_eta, &b_gr_b2_p4_eta);
  ntu->SetBranchAddress("gr_b2_p4_phi", &gr_b2_p4_phi, &b_gr_b2_p4_phi);
  ntu->SetBranchAddress("gr_b2_p4_mass", &gr_b2_p4_mass, &b_gr_b2_p4_mass);
  ntu->SetBranchAddress("gr_hbb_p4_pt", &gr_hbb_p4_pt, &b_gr_hbb_p4_pt);
  ntu->SetBranchAddress("gr_hbb_p4_eta", &gr_hbb_p4_eta, &b_gr_hbb_p4_eta);
  ntu->SetBranchAddress("gr_hbb_p4_phi", &gr_hbb_p4_phi, &b_gr_hbb_p4_phi);
  ntu->SetBranchAddress("gr_hbb_p4_mass", &gr_hbb_p4_mass, &b_gr_hbb_p4_mass);
   
  TLorentzVector* genJet1P4 = new TLorentzVector;
  TLorentzVector* genJet2P4 = new TLorentzVector;
  TLorentzVector* pho1P4 = new TLorentzVector;
  TLorentzVector* pho2P4 = new TLorentzVector;
  TLorentzVector* SumPho = new TLorentzVector();
  TLorentzVector SumJet;
  TLorentzVector Hbb;
 
  std::map<int,TLorentzVector> JetP4;
  
  float nCouples = 0;
  float nCouples_MaxptJets = 0;
  float nCouples_MaxptCouple = 0;
  float nCouples_MaxptOverMCouple = 0;
  float nCouples_MinDeltaMCouple = 0;

  for(int ientry = 0; ientry < ntu->GetEntries(); ientry++){
      if(ientry%1000==0) std::cout<<"--- Reading entry = "<< ientry <<std::endl;
      ntu->GetEntry(ientry);

      int ijet = -1;
      
      genJet1P4->SetPtEtaPhiM(gr_b1_p4_pt,gr_b1_p4_eta,gr_b1_p4_phi,gr_b1_p4_mass);
      genJet2P4->SetPtEtaPhiM(gr_b2_p4_pt,gr_b2_p4_eta,gr_b2_p4_phi,gr_b2_p4_mass);

      Hbb.SetPtEtaPhiM(gr_hbb_p4_pt,gr_hbb_p4_eta,gr_hbb_p4_phi,gr_hbb_p4_mass);
    
      pho1P4->SetPtEtaPhiE(ph1_pt,ph1_eta,ph1_phi,ph1_e);
      pho2P4->SetPtEtaPhiE(ph2_pt,ph2_eta,ph2_phi,ph2_e);
      *SumPho = *pho1P4+*pho2P4;
      float mgg = SumPho->M();
     
      if(isBtagged == 0){
         j1_csvBtag = 1.;
         j2_csvBtag = 1.;
         j3_csvBtag = 1.;
         j4_csvBtag = 1.;
         j5_csvBtag = 1.;
         j6_csvBtag = 1.;
         j7_csvBtag = 1.;
         j8_csvBtag = 1.;
         j9_csvBtag = 1.;
         j10_csvBtag = 1.;
         j11_csvBtag = 1.;
         j12_csvBtag = 1.;
         j13_csvBtag = 1.; 
         j14_csvBtag = 1.; 
         j15_csvBtag = 1.; 
      }
       
      bool passedJet1,passedJet2,passedJet3,passedJet4,passedJet5,passedJet6,passedJet7,passedJet8,passedJet9,passedJet10,
           passedJet11,passedJet12,passedJet13,passedJet14,passedJet15;

      if(SelectJet == 0){
         passedJet1  = true;
         passedJet2  = true;
         passedJet3  = true;
         passedJet4  = true;
         passedJet5  = true;
         passedJet6  = true;
         passedJet7  = true;
         passedJet8  = true;
         passedJet9  = true;
         passedJet10 = true;
         passedJet11 = true;
         passedJet12 = true;
         passedJet13 = true;
         passedJet14 = true;
         passedJet15 = true;
      }else{
         passedJet1  = passCutBasedJetId(j1_pt,j1_eta,j1_betaStarClassic,j1_dR2Mean,nvtx);
         passedJet2  = passCutBasedJetId(j2_pt,j2_eta,j2_betaStarClassic,j2_dR2Mean,nvtx);
         passedJet3  = passCutBasedJetId(j3_pt,j3_eta,j3_betaStarClassic,j3_dR2Mean,nvtx);
         passedJet4  = passCutBasedJetId(j4_pt,j4_eta,j4_betaStarClassic,j4_dR2Mean,nvtx);
         passedJet5  = passCutBasedJetId(j5_pt,j5_eta,j5_betaStarClassic,j5_dR2Mean,nvtx);
         passedJet6  = passCutBasedJetId(j6_pt,j6_eta,j6_betaStarClassic,j6_dR2Mean,nvtx);
         passedJet7  = passCutBasedJetId(j7_pt,j7_eta,j7_betaStarClassic,j7_dR2Mean,nvtx);
         passedJet8  = passCutBasedJetId(j8_pt,j8_eta,j8_betaStarClassic,j8_dR2Mean,nvtx);
         passedJet9  = passCutBasedJetId(j9_pt,j9_eta,j9_betaStarClassic,j9_dR2Mean,nvtx);
         passedJet10 = passCutBasedJetId(j10_pt,j10_eta,j10_betaStarClassic,j10_dR2Mean,nvtx);
         passedJet11 = passCutBasedJetId(j11_pt,j11_eta,j11_betaStarClassic,j11_dR2Mean,nvtx);
         passedJet12 = passCutBasedJetId(j12_pt,j12_eta,j12_betaStarClassic,j12_dR2Mean,nvtx);
         passedJet13 = passCutBasedJetId(j13_pt,j13_eta,j13_betaStarClassic,j13_dR2Mean,nvtx);
         passedJet14 = passCutBasedJetId(j14_pt,j14_eta,j14_betaStarClassic,j14_dR2Mean,nvtx);
         passedJet15 = passCutBasedJetId(j15_pt,j15_eta,j15_betaStarClassic,j15_dR2Mean,nvtx);
      }

      if(j1_e != -1001. && j1_csvBtag > 0.679 && passedJet1 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j1_pt,j1_eta,j1_phi,j1_e); 
      }
 
      if(j2_e != -1001. && j2_csvBtag > 0.679 && passedJet2 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j2_pt,j2_eta,j2_phi,j2_e); 
      }
        
      if(j3_e != -1001. && j3_csvBtag > 0.679 && passedJet3 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j3_pt,j3_eta,j3_phi,j3_e); 
      }

      if(j4_e != -1001. && j4_csvBtag > 0.679 && passedJet4 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j4_pt,j4_eta,j4_phi,j4_e); 
      }

      if(j5_e != -1001. && j5_csvBtag > 0.679 && passedJet5 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j5_pt,j5_eta,j5_phi,j5_e); 
      }

      if(j6_e != -1001. && j6_csvBtag > 0.679 && passedJet6 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j6_pt,j6_eta,j6_phi,j6_e); 
      }   

      if(j7_e != -1001. && j7_csvBtag > 0.679 && passedJet7 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j7_pt,j7_eta,j7_phi,j7_e); 
      }

      if(j8_e != -1001. && j8_csvBtag > 0.679 && passedJet8 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j8_pt,j8_eta,j8_phi,j8_e); 
      }

      if(j9_e != -1001. && j9_csvBtag > 0.679 && passedJet9 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j9_pt,j9_eta,j9_phi,j9_e); 
      }

      if(j10_e != -1001. && j10_csvBtag > 0.679 && passedJet10 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j10_pt,j10_eta,j10_phi,j10_e); 
      }

      if(j11_e != -1001. && j11_csvBtag > 0.679 && passedJet11 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j11_pt,j11_eta,j11_phi,j11_e); 
      }

      if(j12_e != -1001. && j12_csvBtag > 0.679 && passedJet12 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j12_pt,j12_eta,j12_phi,j12_e); 
      }

      if(j13_e != -1001. && j13_csvBtag > 0.679 && passedJet13 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j13_pt,j13_eta,j13_phi,j13_e); 
      }

      if(j14_e != -1001. && j14_csvBtag > 0.679 && passedJet14 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j14_pt,j14_eta,j14_phi,j14_e); 
      }

      if(j15_e != -1001. && j5_csvBtag > 0.679 && passedJet15 == true){
         ijet++;
         JetP4[ijet].SetPtEtaPhiE(j15_pt,j15_eta,j15_phi,j15_e); 
      }

      if(JetP4.size() < 2) continue;
      if(evweight == 0.) continue;

      nCouples = nCouples+evweight;

      //MaxptJets
      int ijet1 = chooseMaxptJets(&JetP4).at(0);
      int ijet2 = chooseMaxptJets(&JetP4).at(1);
      SumJet = JetP4.at(ijet1) + JetP4.at(ijet2);
      if((JetP4.at(ijet1).DeltaR(*genJet1P4) < dR && JetP4.at(ijet2).DeltaR(*genJet2P4) < dR) || (JetP4.at(ijet2).DeltaR(*genJet1P4) < dR && JetP4.at(ijet1).DeltaR(*genJet2P4) < dR)) nCouples_MaxptJets =nCouples_MaxptJets+evweight;
      //if(Hbb.DeltaR(SumJet) < dR) nCouples_MaxptJets =nCouples_MaxptJets+evweight;
     
      //MaxptCouple
      ijet1 = chooseMaxptCouple(&JetP4).at(0);
      ijet2 = chooseMaxptCouple(&JetP4).at(1);
      SumJet = JetP4.at(ijet1) + JetP4.at(ijet2);
      if((JetP4.at(ijet1).DeltaR(*genJet1P4) < dR && JetP4.at(ijet2).DeltaR(*genJet2P4) < dR) || (JetP4.at(ijet2).DeltaR(*genJet1P4) < dR && JetP4.at(ijet1).DeltaR(*genJet2P4) < dR)) nCouples_MaxptCouple = nCouples_MaxptCouple+evweight;
      //if(Hbb.DeltaR(SumJet) < dR) nCouples_MaxptCouple = nCouples_MaxptCouple+evweight;

      //MaxptOverMCouple
      ijet1 = chooseMaxptOverMCouple(&JetP4).at(0);
      ijet2 = chooseMaxptOverMCouple(&JetP4).at(1);
      SumJet = JetP4.at(ijet1) + JetP4.at(ijet2);
      if((JetP4.at(ijet1).DeltaR(*genJet1P4) < dR && JetP4.at(ijet2).DeltaR(*genJet2P4) < dR) || (JetP4.at(ijet2).DeltaR(*genJet1P4) < dR && JetP4.at(ijet1).DeltaR(*genJet2P4) < dR)) nCouples_MaxptOverMCouple = nCouples_MaxptOverMCouple+evweight;
      //if(Hbb.DeltaR(SumJet) < dR) nCouples_MaxptOverMCouple = nCouples_MaxptOverMCouple+evweight;

      //MinDeltaMCouple
      ijet1 = chooseMinDeltaMCouple(&JetP4,mgg).at(0);
      ijet2 = chooseMinDeltaMCouple(&JetP4,mgg).at(1); 
      SumJet = JetP4.at(ijet1) + JetP4.at(ijet2);
      if((JetP4.at(ijet1).DeltaR(*genJet1P4) < dR && JetP4.at(ijet2).DeltaR(*genJet2P4) < dR) || (JetP4.at(ijet2).DeltaR(*genJet1P4) < dR && JetP4.at(ijet1).DeltaR(*genJet2P4) < dR)) nCouples_MinDeltaMCouple = nCouples_MinDeltaMCouple+evweight;
      //if(Hbb.DeltaR(SumJet) < dR) nCouples_MinDeltaMCouple = nCouples_MinDeltaMCouple+evweight;
   
      JetP4.clear();

  }

  std::cout << "Eff_MaxptJets = " << float(nCouples_MaxptJets)/float(nCouples)*100. << " %" << std::endl;
  std::cout << "Eff_MaxptCouple = " << float(nCouples_MaxptCouple)/float(nCouples)*100. << " %" << std::endl;
  std::cout << "Eff_MaxptOverMCouple = " << float(nCouples_MaxptOverMCouple)/float(nCouples)*100. << " %" << std::endl;
  std::cout << "Eff_MinDeltaMCouple = " << float(nCouples_MinDeltaMCouple)/float(nCouples)*100. << " %" << std::endl;
}

bool passCutBasedJetId(const float& pt,const float& eta, const float& betaStarClassic, const float& dR2Mean, const int& nvtx){

   bool isGood = false;
  
   if(pt>25. && fabs(eta)<2.5 && betaStarClassic<=0.2*log(nvtx-0.64) && dR2Mean<=0.06) isGood = true;

   return isGood;  
}

std::vector<int> chooseMaxptJets(std::map<int,TLorentzVector>* JetP4){

  int jet1 = -1;
  int jet2 = -1; 
   
  float maxPt = 0.;
  
  for(unsigned int ii = 0; ii <  JetP4->size(); ii++){

      if(JetP4->at(ii).Pt() > maxPt){
             maxPt = JetP4->at(ii).Pt();
             jet1 = ii;
      }       
  }
  
  maxPt = 0.;
  
  for(unsigned int jj = 0; jj <  JetP4->size(); jj++){

      if(int(jj) == int(jet1)) continue;

      if(JetP4->at(jj).Pt() > maxPt){

             maxPt = JetP4->at(jj).Pt();
             jet2 = jj;
      }       
  }

  std::vector<int> outvec;
  outvec.push_back(jet1);
  outvec.push_back(jet2);

  return outvec;

}


std::vector<int> chooseMaxptCouple(std::map<int,TLorentzVector>* JetP4){

  int jet1 = -1;
  int jet2 = -1; 
   
  float maxPt = 0.;
  TLorentzVector Sum;
  
  for(unsigned int ii = 0; ii <  JetP4->size(); ii++){

      for(unsigned int jj = 0; jj <  JetP4->size(); jj++){

          if(ii == jj) continue;  
          
          Sum = JetP4->at(ii)+JetP4->at(jj);
          
          if(Sum.Pt() > maxPt){

             maxPt = Sum.Pt();
             jet1 = ii;
             jet2 = jj;

          }     
      }  
  }

  std::vector<int> outvec;
  outvec.push_back(jet1);
  outvec.push_back(jet2);

  return outvec;

}

std::vector<int> chooseMaxptOverMCouple(std::map<int,TLorentzVector>* JetP4){

  int jet1 = -1;
  int jet2 = -1; 
   
  float maxPtOverM = 0.;
  TLorentzVector Sum;
  
  for(unsigned int ii = 0; ii <  JetP4->size(); ii++){

      for(unsigned int jj = 0; jj <  JetP4->size(); jj++){

          if(ii == jj) continue;  
          
          Sum = JetP4->at(ii)+JetP4->at(ii);
          
          if(Sum.Pt()/Sum.M() > maxPtOverM){

             maxPtOverM = Sum.Pt()/Sum.M();
             jet1 = ii;
             jet2 = jj;

          }     
      }  
  }

  std::vector<int> outvec;
  outvec.push_back(jet1);
  outvec.push_back(jet2);

  return outvec;

}

std::vector<int> chooseMinDeltaMCouple(std::map<int,TLorentzVector>* JetP4,const float mgg){

  int jet1 = -1;
  int jet2 = -1; 
   
  float delltaM = 0.;
  TLorentzVector Sum;
  
  for(unsigned int ii = 0; ii <  JetP4->size(); ii++){

      for(unsigned int jj = 0; jj <  JetP4->size(); jj++){

          if(ii == jj) continue;  
          
          Sum = JetP4->at(ii)+JetP4->at(jj);
          
          if(1./fabs(Sum.M()-mgg) > delltaM){

             delltaM = 1./fabs(Sum.M()-mgg);
             jet1 = ii;
             jet2 = jj;
          }     
      }  
  }

  std::vector<int> outvec;
  outvec.push_back(jet1);
  outvec.push_back(jet2);

  return outvec;

}



