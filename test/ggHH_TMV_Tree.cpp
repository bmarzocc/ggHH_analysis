
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"
#include "BTagUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
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
#include "THStack.h"
#include "TBranch.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

int main(int argc, char** argv)
{
   
  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> Analyses::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");
  std::string inputRef = gConfigParser -> readStringOption("Input::inputRef");
  std::string inputTreeRef = gConfigParser -> readStringOption("Input::inputTreeRef");
  std::string inputData = gConfigParser -> readStringOption("Input::inputData");
  std::string inputTreeData = gConfigParser -> readStringOption("Input::inputTreeData");

  std::string outputDir = gConfigParser -> readStringOption("Output::outputDir");

  TChain* ntu_ref = new TChain(inputTreeRef.c_str());
  ntu_ref->Add(inputRef.c_str());

  TChain* ntu_data = new TChain(inputTreeData.c_str());
  ntu_data->Add(inputData.c_str());

  std::map<int,TChain*> ntu;
  std::map<int, TFile*> Files;
   
  char fileName[2000];
  char trees[500];
  FILE *f_trees;
  f_trees = fopen(inputList.c_str(),"r");
  int pos = 0;
  int pos_total = 0;
  
  while(fscanf(f_trees,"%s %s\n", fileName, trees) !=EOF ){

    std::string FILENAME = std::string(fileName);
    if(FILENAME.find("#") != std::string::npos) continue;
    std::cout << "\nReading input: " << fileName << " - " << trees << std::endl;
    ntu[pos] = new TChain(trees);
    ntu[pos] -> Add(fileName);
   
    pos++;
    
  }
  
  pos_total = pos;
  
  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file: " << ntu[ii]->GetName() << " is empty" << std::endl; 
         return -1;
      }
  }

  int           category;
  float         event;
  float         evweight;
  float         pho1_pt;
  float         pho1_e;
  float         pho1_phi;
  float         pho1_eta;
  float         pho1_mass;
  float         pho1_r9;
  float         pho1_sieie;
  float         pho1_hoe;
  float         pho2_pt;
  float         pho2_e;
  float         pho2_phi;
  float         pho2_eta;
  float         pho2_mass;
  float         pho2_r9;
  float         pho2_sieie;
  float         pho2_hoe;
  float         jet1_pt;
  float         jet1_e;
  float         jet1_phi;
  float         jet1_eta;
  float         jet1_mass;
  float         jet1_csvBtag;
  float         jet1_btagSF_M;
  float         jet1_btagSFErrorUp_M;
  float         jet1_btagSFErrorDown_M;
  float         jet1_btagEff_M;
  float         jet1_btagEffError_M;
  int           jet1_flavour;
  float         jet1_betaStarClassic;
  float         jet1_dR2Mean;
  float         jet2_pt;
  float         jet2_e;
  float         jet2_phi;
  float         jet2_eta;
  float         jet2_mass;
  float         jet2_csvBtag;
  float         jet2_btagSF_M;
  float         jet2_btagSFErrorUp_M;
  float         jet2_btagSFErrorDown_M;
  float         jet2_btagEff_M;
  float         jet2_btagEffError_M;
  int           jet2_flavour;
  float         jet2_betaStarClassic;
  float         jet2_dR2Mean;
  float         jj_pt;
  float         jj_e;
  float         jj_phi;
  float         jj_eta;
  float         jj_mass;
  float         jj_DR;
  float         gg_pt;
  float         gg_e;
  float         gg_phi;
  float         gg_eta;
  float         gg_mass;
  float         costhetastar;
  float         minDRgj;
  float         ggjj_pt;
  float         ggjj_phi;
  float         ggjj_eta;
  float         ggjj_mass;

  // List of branches
   TBranch        *b_category;   //!
   TBranch        *b_event;   //!
   TBranch        *b_evweight;   //!
   TBranch        *b_pho1_pt;   //!
   TBranch        *b_pho1_e;   //!
   TBranch        *b_pho1_phi;   //!
   TBranch        *b_pho1_eta;   //!
   TBranch        *b_pho1_mass;   //!
   TBranch        *b_pho1_r9;   //!
   TBranch        *b_pho1_sieie;   //!
   TBranch        *b_pho1_hoe;   //!
   TBranch        *b_pho2_pt;   //!
   TBranch        *b_pho2_e;   //!
   TBranch        *b_pho2_phi;   //!
   TBranch        *b_pho2_eta;   //!
   TBranch        *b_pho2_mass;   //!
   TBranch        *b_pho2_r9;   //!
   TBranch        *b_pho2_sieie;   //!
   TBranch        *b_pho2_hoe;   //!
   TBranch        *b_jet1_pt;   //!
   TBranch        *b_jet1_e;   //!
   TBranch        *b_jet1_phi;   //!
   TBranch        *b_jet1_eta;   //!
   TBranch        *b_jet1_mass;   //!
   TBranch        *b_jet1_csvBtag;   //!
   TBranch        *b_jet1_btagSF_M;   //!
   TBranch        *b_jet1_btagSFErrorUp_M;   //!
   TBranch        *b_jet1_btagSFErrorDown_M;   //!
   TBranch        *b_jet1_btagEff_M;   //!
   TBranch        *b_jet1_btagEffError_M;   //!
   TBranch        *b_jet1_flavour;   //!
   TBranch        *b_jet1_betaStarClassic;   //!
   TBranch        *b_jet1_dR2Mean;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet2_e;   //!
   TBranch        *b_jet2_phi;   //!
   TBranch        *b_jet2_eta;   //!
   TBranch        *b_jet2_mass;   //!
   TBranch        *b_jet2_csvBtag;   //!
   TBranch        *b_jet2_btagSF_M;   //!
   TBranch        *b_jet2_btagSFErrorUp_M;   //!
   TBranch        *b_jet2_btagSFErrorDown_M;   //!
   TBranch        *b_jet2_btagEff_M;   //!
   TBranch        *b_jet2_btagEffError_M;   //!
   TBranch        *b_jet2_flavour;   //!
   TBranch        *b_jet2_betaStarClassic;   //!
   TBranch        *b_jet2_dR2Mean;   //!
   TBranch        *b_jj_pt;   //!
   TBranch        *b_jj_e;   //!
   TBranch        *b_jj_phi;   //!
   TBranch        *b_jj_eta;   //!
   TBranch        *b_jj_mass;   //!
   TBranch        *b_jj_DR;   //!
   TBranch        *b_gg_pt;   //!
   TBranch        *b_gg_e;   //!
   TBranch        *b_gg_phi;   //!
   TBranch        *b_gg_eta;   //!
   TBranch        *b_gg_mass;   //!
   TBranch        *b_costhetastar;   //!
   TBranch        *b_minDRgj;   //!
   TBranch        *b_ggjj_pt;   //!
   TBranch        *b_ggjj_phi;   //!
   TBranch        *b_ggjj_eta;   //!
   TBranch        *b_ggjj_mass;   //!

  for(int ii = 0; ii < pos_total; ii++){
   ntu[ii]->SetBranchAddress("category", &category, &b_category);
   ntu[ii]->SetBranchAddress("event", &event, &b_event);
   ntu[ii]->SetBranchAddress("evweight", &evweight, &b_evweight);
   ntu[ii]->SetBranchAddress("pho1_pt", &pho1_pt, &b_pho1_pt);
   ntu[ii]->SetBranchAddress("pho1_e", &pho1_e, &b_pho1_e);
   ntu[ii]->SetBranchAddress("pho1_phi", &pho1_phi, &b_pho1_phi);
   ntu[ii]->SetBranchAddress("pho1_eta", &pho1_eta, &b_pho1_eta);
   ntu[ii]->SetBranchAddress("pho1_mass", &pho1_mass, &b_pho1_mass);
   ntu[ii]->SetBranchAddress("pho1_r9", &pho1_r9, &b_pho1_r9);
   ntu[ii]->SetBranchAddress("pho1_sieie", &pho1_sieie, &b_pho1_sieie);
   ntu[ii]->SetBranchAddress("pho1_hoe", &pho1_hoe, &b_pho1_hoe);
   ntu[ii]->SetBranchAddress("pho2_pt", &pho2_pt, &b_pho2_pt);
   ntu[ii]->SetBranchAddress("pho2_e", &pho2_e, &b_pho2_e);
   ntu[ii]->SetBranchAddress("pho2_phi", &pho2_phi, &b_pho2_phi);
   ntu[ii]->SetBranchAddress("pho2_eta", &pho2_eta, &b_pho2_eta);
   ntu[ii]->SetBranchAddress("pho2_mass", &pho2_mass, &b_pho2_mass);
   ntu[ii]->SetBranchAddress("pho2_r9", &pho2_r9, &b_pho2_r9);
   ntu[ii]->SetBranchAddress("pho2_sieie", &pho2_sieie, &b_pho2_sieie);
   ntu[ii]->SetBranchAddress("pho2_hoe", &pho2_hoe, &b_pho2_hoe);
   ntu[ii]->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   ntu[ii]->SetBranchAddress("jet1_e", &jet1_e, &b_jet1_e);
   ntu[ii]->SetBranchAddress("jet1_phi", &jet1_phi, &b_jet1_phi);
   ntu[ii]->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   ntu[ii]->SetBranchAddress("jet1_mass", &jet1_mass, &b_jet1_mass);
   ntu[ii]->SetBranchAddress("jet1_csvBtag", &jet1_csvBtag, &b_jet1_csvBtag);
   ntu[ii]->SetBranchAddress("jet1_btagSF_M", &jet1_btagSF_M, &b_jet1_btagSF_M);
   ntu[ii]->SetBranchAddress("jet1_btagSFErrorUp_M", &jet1_btagSFErrorUp_M, &b_jet1_btagSFErrorUp_M);
   ntu[ii]->SetBranchAddress("jet1_btagSFErrorDown_M", &jet1_btagSFErrorDown_M, &b_jet1_btagSFErrorDown_M);
   ntu[ii]->SetBranchAddress("jet1_btagEff_M", &jet1_btagEff_M, &b_jet1_btagEff_M);
   ntu[ii]->SetBranchAddress("jet1_btagEffError_M", &jet1_btagEffError_M, &b_jet1_btagEffError_M);
   ntu[ii]->SetBranchAddress("jet1_flavour", &jet1_flavour, &b_jet1_flavour);
   ntu[ii]->SetBranchAddress("jet1_betaStarClassic", &jet1_betaStarClassic, &b_jet1_betaStarClassic);
   ntu[ii]->SetBranchAddress("jet1_dR2Mean", &jet1_dR2Mean, &b_jet1_dR2Mean);
   ntu[ii]->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   ntu[ii]->SetBranchAddress("jet2_e", &jet2_e, &b_jet2_e);
   ntu[ii]->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   ntu[ii]->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   ntu[ii]->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   ntu[ii]->SetBranchAddress("jet2_csvBtag", &jet2_csvBtag, &b_jet2_csvBtag);
   ntu[ii]->SetBranchAddress("jet2_btagSF_M", &jet2_btagSF_M, &b_jet2_btagSF_M);
   ntu[ii]->SetBranchAddress("jet2_btagSFErrorUp_M", &jet2_btagSFErrorUp_M, &b_jet2_btagSFErrorUp_M);
   ntu[ii]->SetBranchAddress("jet2_btagSFErrorDown_M", &jet2_btagSFErrorDown_M, &b_jet2_btagSFErrorDown_M);
   ntu[ii]->SetBranchAddress("jet2_btagEff_M", &jet2_btagEff_M, &b_jet2_btagEff_M);
   ntu[ii]->SetBranchAddress("jet2_btagEffError_M", &jet2_btagEffError_M, &b_jet2_btagEffError_M);
   ntu[ii]->SetBranchAddress("jet2_flavour", &jet2_flavour, &b_jet2_flavour);
   ntu[ii]->SetBranchAddress("jet2_betaStarClassic", &jet2_betaStarClassic, &b_jet2_betaStarClassic);
   ntu[ii]->SetBranchAddress("jet2_dR2Mean", &jet2_dR2Mean, &b_jet2_dR2Mean);
   ntu[ii]->SetBranchAddress("jj_pt", &jj_pt, &b_jj_pt);
   ntu[ii]->SetBranchAddress("jj_e", &jj_e, &b_jj_e);
   ntu[ii]->SetBranchAddress("jj_phi", &jj_phi, &b_jj_phi);
   ntu[ii]->SetBranchAddress("jj_eta", &jj_eta, &b_jj_eta);
   ntu[ii]->SetBranchAddress("jj_mass", &jj_mass, &b_jj_mass);
   ntu[ii]->SetBranchAddress("jj_DR", &jj_DR, &b_jj_DR);
   ntu[ii]->SetBranchAddress("gg_pt", &gg_pt, &b_gg_pt);
   ntu[ii]->SetBranchAddress("gg_e", &gg_e, &b_gg_e);
   ntu[ii]->SetBranchAddress("gg_phi", &gg_phi, &b_gg_phi);
   ntu[ii]->SetBranchAddress("gg_eta", &gg_eta, &b_gg_eta);
   ntu[ii]->SetBranchAddress("gg_mass", &gg_mass, &b_gg_mass);
   ntu[ii]->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   ntu[ii]->SetBranchAddress("minDRgj", &minDRgj, &b_minDRgj);
   ntu[ii]->SetBranchAddress("ggjj_pt", &ggjj_pt, &b_ggjj_pt);
   ntu[ii]->SetBranchAddress("ggjj_eta", &ggjj_eta, &b_ggjj_eta);
   ntu[ii]->SetBranchAddress("ggjj_phi", &ggjj_phi, &b_ggjj_phi);
   ntu[ii]->SetBranchAddress("ggjj_mass", &ggjj_mass, &b_ggjj_mass);
  }

   ntu_ref->SetBranchAddress("category", &category, &b_category);
   ntu_ref->SetBranchAddress("event", &event, &b_event);
   ntu_ref->SetBranchAddress("evweight", &evweight, &b_evweight);
   ntu_ref->SetBranchAddress("pho1_pt", &pho1_pt, &b_pho1_pt);
   ntu_ref->SetBranchAddress("pho1_e", &pho1_e, &b_pho1_e);
   ntu_ref->SetBranchAddress("pho1_phi", &pho1_phi, &b_pho1_phi);
   ntu_ref->SetBranchAddress("pho1_eta", &pho1_eta, &b_pho1_eta);
   ntu_ref->SetBranchAddress("pho1_mass", &pho1_mass, &b_pho1_mass);
   ntu_ref->SetBranchAddress("pho1_r9", &pho1_r9, &b_pho1_r9);
   ntu_ref->SetBranchAddress("pho1_sieie", &pho1_sieie, &b_pho1_sieie);
   ntu_ref->SetBranchAddress("pho1_hoe", &pho1_hoe, &b_pho1_hoe);
   ntu_ref->SetBranchAddress("pho2_pt", &pho2_pt, &b_pho2_pt);
   ntu_ref->SetBranchAddress("pho2_e", &pho2_e, &b_pho2_e);
   ntu_ref->SetBranchAddress("pho2_phi", &pho2_phi, &b_pho2_phi);
   ntu_ref->SetBranchAddress("pho2_eta", &pho2_eta, &b_pho2_eta);
   ntu_ref->SetBranchAddress("pho2_mass", &pho2_mass, &b_pho2_mass);
   ntu_ref->SetBranchAddress("pho2_r9", &pho2_r9, &b_pho2_r9);
   ntu_ref->SetBranchAddress("pho2_sieie", &pho2_sieie, &b_pho2_sieie);
   ntu_ref->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   ntu_ref->SetBranchAddress("jet1_e", &jet1_e, &b_jet1_e);
   ntu_ref->SetBranchAddress("jet1_phi", &jet1_phi, &b_jet1_phi);
   ntu_ref->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   ntu_ref->SetBranchAddress("jet1_mass", &jet1_mass, &b_jet1_mass);
   ntu_ref->SetBranchAddress("jet1_csvBtag", &jet1_csvBtag, &b_jet1_csvBtag);
   ntu_ref->SetBranchAddress("jet1_btagSF_M", &jet1_btagSF_M, &b_jet1_btagSF_M);
   ntu_ref->SetBranchAddress("jet1_btagSFErrorUp_M", &jet1_btagSFErrorUp_M, &b_jet1_btagSFErrorUp_M);
   ntu_ref->SetBranchAddress("jet1_btagSFErrorDown_M", &jet1_btagSFErrorDown_M, &b_jet1_btagSFErrorDown_M);
   ntu_ref->SetBranchAddress("jet1_btagEff_M", &jet1_btagEff_M, &b_jet1_btagEff_M);
   ntu_ref->SetBranchAddress("jet1_btagEffError_M", &jet1_btagEffError_M, &b_jet1_btagEffError_M);
   ntu_ref->SetBranchAddress("jet1_flavour", &jet1_flavour, &b_jet1_flavour);
   ntu_ref->SetBranchAddress("jet1_betaStarClassic", &jet1_betaStarClassic, &b_jet1_betaStarClassic);
   ntu_ref->SetBranchAddress("jet1_dR2Mean", &jet1_dR2Mean, &b_jet1_dR2Mean);
   ntu_ref->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   ntu_ref->SetBranchAddress("jet2_e", &jet2_e, &b_jet2_e);
   ntu_ref->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   ntu_ref->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   ntu_ref->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   ntu_ref->SetBranchAddress("jet2_csvBtag", &jet2_csvBtag, &b_jet2_csvBtag);
   ntu_ref->SetBranchAddress("jet2_btagSF_M", &jet2_btagSF_M, &b_jet2_btagSF_M);
   ntu_ref->SetBranchAddress("jet2_btagSFErrorUp_M", &jet2_btagSFErrorUp_M, &b_jet2_btagSFErrorUp_M);
   ntu_ref->SetBranchAddress("jet2_btagSFErrorDown_M", &jet2_btagSFErrorDown_M, &b_jet2_btagSFErrorDown_M);
   ntu_ref->SetBranchAddress("jet2_btagEff_M", &jet2_btagEff_M, &b_jet2_btagEff_M);
   ntu_ref->SetBranchAddress("jet2_btagEffError_M", &jet2_btagEffError_M, &b_jet2_btagEffError_M);
   ntu_ref->SetBranchAddress("jet2_flavour", &jet2_flavour, &b_jet2_flavour);
   ntu_ref->SetBranchAddress("jet2_betaStarClassic", &jet2_betaStarClassic, &b_jet2_betaStarClassic);
   ntu_ref->SetBranchAddress("jet2_dR2Mean", &jet2_dR2Mean, &b_jet2_dR2Mean);
   ntu_ref->SetBranchAddress("jj_pt", &jj_pt, &b_jj_pt);
   ntu_ref->SetBranchAddress("jj_e", &jj_e, &b_jj_e);
   ntu_ref->SetBranchAddress("jj_phi", &jj_phi, &b_jj_phi);
   ntu_ref->SetBranchAddress("jj_eta", &jj_eta, &b_jj_eta);
   ntu_ref->SetBranchAddress("jj_mass", &jj_mass, &b_jj_mass);
   ntu_ref->SetBranchAddress("jj_DR", &jj_DR, &b_jj_DR);
   ntu_ref->SetBranchAddress("gg_pt", &gg_pt, &b_gg_pt);
   ntu_ref->SetBranchAddress("gg_e", &gg_e, &b_gg_e);
   ntu_ref->SetBranchAddress("gg_phi", &gg_phi, &b_gg_phi);
   ntu_ref->SetBranchAddress("gg_eta", &gg_eta, &b_gg_eta);
   ntu_ref->SetBranchAddress("gg_mass", &gg_mass, &b_gg_mass);
   ntu_ref->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   ntu_ref->SetBranchAddress("minDRgj", &minDRgj, &b_minDRgj);
   ntu_ref->SetBranchAddress("ggjj_pt", &ggjj_pt, &b_ggjj_pt);
   ntu_ref->SetBranchAddress("ggjj_eta", &ggjj_eta, &b_ggjj_eta);
   ntu_ref->SetBranchAddress("ggjj_phi", &ggjj_phi, &b_ggjj_phi);
   ntu_ref->SetBranchAddress("ggjj_mass", &ggjj_mass, &b_ggjj_mass);

   ntu_data->SetBranchAddress("event", &event, &b_event);
   ntu_data->SetBranchAddress("category", &category, &b_category);
   ntu_data->SetBranchAddress("evweight", &evweight, &b_evweight);
   ntu_data->SetBranchAddress("pho1_pt", &pho1_pt, &b_pho1_pt);
   ntu_data->SetBranchAddress("pho1_e", &pho1_e, &b_pho1_e);
   ntu_data->SetBranchAddress("pho1_phi", &pho1_phi, &b_pho1_phi);
   ntu_data->SetBranchAddress("pho1_eta", &pho1_eta, &b_pho1_eta);
   ntu_data->SetBranchAddress("pho1_mass", &pho1_mass, &b_pho1_mass);
   ntu_data->SetBranchAddress("pho1_r9", &pho1_r9, &b_pho1_r9);
   ntu_data->SetBranchAddress("pho1_sieie", &pho1_sieie, &b_pho1_sieie);
   ntu_data->SetBranchAddress("pho1_hoe", &pho1_hoe, &b_pho1_hoe);
   ntu_data->SetBranchAddress("pho2_pt", &pho2_pt, &b_pho2_pt);
   ntu_data->SetBranchAddress("pho2_e", &pho2_e, &b_pho2_e);
   ntu_data->SetBranchAddress("pho2_phi", &pho2_phi, &b_pho2_phi);
   ntu_data->SetBranchAddress("pho2_eta", &pho2_eta, &b_pho2_eta);
   ntu_data->SetBranchAddress("pho2_mass", &pho2_mass, &b_pho2_mass);
   ntu_data->SetBranchAddress("pho2_r9", &pho2_r9, &b_pho2_r9);
   ntu_data->SetBranchAddress("pho2_sieie", &pho2_sieie, &b_pho2_sieie);
   ntu_data->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   ntu_data->SetBranchAddress("jet1_e", &jet1_e, &b_jet1_e);
   ntu_data->SetBranchAddress("jet1_phi", &jet1_phi, &b_jet1_phi);
   ntu_data->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   ntu_data->SetBranchAddress("jet1_mass", &jet1_mass, &b_jet1_mass);
   ntu_data->SetBranchAddress("jet1_csvBtag", &jet1_csvBtag, &b_jet1_csvBtag);
   ntu_data->SetBranchAddress("jet1_btagSF_M", &jet1_btagSF_M, &b_jet1_btagSF_M);
   ntu_data->SetBranchAddress("jet1_btagSFErrorUp_M", &jet1_btagSFErrorUp_M, &b_jet1_btagSFErrorUp_M);
   ntu_data->SetBranchAddress("jet1_btagSFErrorDown_M", &jet1_btagSFErrorDown_M, &b_jet1_btagSFErrorDown_M);
   ntu_data->SetBranchAddress("jet1_btagEff_M", &jet1_btagEff_M, &b_jet1_btagEff_M);
   ntu_data->SetBranchAddress("jet1_btagEffError_M", &jet1_btagEffError_M, &b_jet1_btagEffError_M);
   ntu_data->SetBranchAddress("jet1_flavour", &jet1_flavour, &b_jet1_flavour);
   ntu_data->SetBranchAddress("jet1_betaStarClassic", &jet1_betaStarClassic, &b_jet1_betaStarClassic);
   ntu_data->SetBranchAddress("jet1_dR2Mean", &jet1_dR2Mean, &b_jet1_dR2Mean);
   ntu_data->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   ntu_data->SetBranchAddress("jet2_e", &jet2_e, &b_jet2_e);
   ntu_data->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   ntu_data->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   ntu_data->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   ntu_data->SetBranchAddress("jet2_csvBtag", &jet2_csvBtag, &b_jet2_csvBtag);
   ntu_data->SetBranchAddress("jet2_btagSF_M", &jet2_btagSF_M, &b_jet2_btagSF_M);
   ntu_data->SetBranchAddress("jet2_btagSFErrorUp_M", &jet2_btagSFErrorUp_M, &b_jet2_btagSFErrorUp_M);
   ntu_data->SetBranchAddress("jet2_btagSFErrorDown_M", &jet2_btagSFErrorDown_M, &b_jet2_btagSFErrorDown_M);
   ntu_data->SetBranchAddress("jet2_btagEff_M", &jet2_btagEff_M, &b_jet2_btagEff_M);
   ntu_data->SetBranchAddress("jet2_btagEffError_M", &jet2_btagEffError_M, &b_jet2_btagEffError_M);
   ntu_data->SetBranchAddress("jet2_flavour", &jet2_flavour, &b_jet2_flavour);
   ntu_data->SetBranchAddress("jet2_betaStarClassic", &jet2_betaStarClassic, &b_jet2_betaStarClassic);
   ntu_data->SetBranchAddress("jet2_dR2Mean", &jet2_dR2Mean, &b_jet2_dR2Mean);
   ntu_data->SetBranchAddress("jj_pt", &jj_pt, &b_jj_pt);
   ntu_data->SetBranchAddress("jj_e", &jj_e, &b_jj_e);
   ntu_data->SetBranchAddress("jj_phi", &jj_phi, &b_jj_phi);
   ntu_data->SetBranchAddress("jj_eta", &jj_eta, &b_jj_eta);
   ntu_data->SetBranchAddress("jj_mass", &jj_mass, &b_jj_mass);
   ntu_data->SetBranchAddress("jj_DR", &jj_DR, &b_jj_DR);
   ntu_data->SetBranchAddress("gg_pt", &gg_pt, &b_gg_pt);
   ntu_data->SetBranchAddress("gg_e", &gg_e, &b_gg_e);
   ntu_data->SetBranchAddress("gg_phi", &gg_phi, &b_gg_phi);
   ntu_data->SetBranchAddress("gg_eta", &gg_eta, &b_gg_eta);
   ntu_data->SetBranchAddress("gg_mass", &gg_mass, &b_gg_mass);
   ntu_data->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   ntu_data->SetBranchAddress("minDRgj", &minDRgj, &b_minDRgj);
   ntu_data->SetBranchAddress("ggjj_pt", &ggjj_pt, &b_ggjj_pt);
   ntu_data->SetBranchAddress("ggjj_eta", &ggjj_eta, &b_ggjj_eta);
   ntu_data->SetBranchAddress("ggjj_phi", &ggjj_phi, &b_ggjj_phi);
   ntu_data->SetBranchAddress("ggjj_mass", &ggjj_mass, &b_ggjj_mass);

  int           category_output;
  float         event_output;
  float         evweight_output;
  float         pho1_pt_output;
  float         pho1_e_output;
  float         pho1_phi_output;
  float         pho1_eta_output;
  float         pho1_mass_output;
  float         pho1_r9_output;
  float         pho1_sieie_output;
  float         pho1_hoe_output;
  float         pho2_pt_output;
  float         pho2_e_output;
  float         pho2_phi_output;
  float         pho2_eta_output;
  float         pho2_mass_output;
  float         pho2_r9_output;
  float         pho2_sieie_output;
  float         pho2_hoe_output;
  float         jet1_pt_output;
  float         jet1_e_output;
  float         jet1_phi_output;
  float         jet1_eta_output;
  float         jet1_mass_output;
  float         jet1_csvBtag_output;
  float         jet1_btagSF_M_output;
  float         jet1_btagSFErrorUp_M_output;
  float         jet1_btagSFErrorDown_M_output;
  float         jet1_btagEff_M_output;
  float         jet1_btagEffError_M_output;
  int           jet1_flavour_output;
  float         jet1_betaStarClassic_output;
  float         jet1_dR2Mean_output;
  float         jet2_pt_output;
  float         jet2_e_output;
  float         jet2_phi_output;
  float         jet2_eta_output;
  float         jet2_mass_output;
  float         jet2_csvBtag_output;
  float         jet2_btagSF_M_output;
  float         jet2_btagSFErrorUp_M_output;
  float         jet2_btagSFErrorDown_M_output;
  float         jet2_btagEff_M_output;
  float         jet2_btagEffError_M_output;
  int           jet2_flavour_output;
  float         jet2_betaStarClassic_output;
  float         jet2_dR2Mean_output;
  float         jj_pt_output;
  float         jj_e_output;
  float         jj_phi_output;
  float         jj_eta_output;
  float         jj_mass_output;
  float         jj_DR_output;
  float         gg_pt_output;
  float         gg_e_output;
  float         gg_phi_output;
  float         gg_eta_output;
  float         gg_mass_output;
  float         gg_DR_output;
  float         costhetastar_output;
  float         minDRgj_output;
  float         ggjj_pt_output;
  float         ggjj_e_output;
  float         ggjj_phi_output;
  float         ggjj_eta_output;
  float         ggjj_mass_output;
  float         ggjj_DR_output;
   
  std::map<int,TTree*> outTree;
  outTree[0] = new TTree("TreeBkg","TreeBkg");
  outTree[1] = new TTree("TreeSignal","TreeSignal");
  outTree[2] = new TTree("TreeData","TreeData");
  
  for(unsigned int ii = 0; ii < 3; ii++){

   outTree[ii]->SetDirectory(0);

   outTree[ii] -> Branch("event",                   &event_output,                       "event/F");
   outTree[ii] -> Branch("category",                &category_output,                    "category/I");
   outTree[ii] -> Branch("evweight",                &evweight_output,                    "evweight/F");
   outTree[ii] -> Branch("pho1_pt",                 &pho1_pt_output,                     "pho1_pt/F");
   outTree[ii] -> Branch("pho1_e",                  &pho1_e_output,                      "pho1_e/F");
   outTree[ii] -> Branch("pho1_eta",                &pho1_eta_output,                    "pho1_eta/F");
   outTree[ii] -> Branch("pho1_phi",                &pho1_phi_output,                    "pho1_phi/F");
   outTree[ii] -> Branch("pho1_mass",               &pho1_mass_output,                   "pho1_mass/F");
   outTree[ii] -> Branch("pho1_sieie",              &pho1_sieie_output,                  "pho1_sieie/F");
   outTree[ii] -> Branch("pho1_hoe",                &pho1_hoe_output,                    "pho1_hoe/F");
   outTree[ii] -> Branch("pho1_r9",                 &pho1_r9_output,                     "pho1_r9/F");
   outTree[ii] -> Branch("pho2_pt",                 &pho2_pt_output,                     "pho2_pt/F");
   outTree[ii] -> Branch("pho2_e",                  &pho2_e_output,                      "pho2_e/F");
   outTree[ii] -> Branch("pho2_eta",                &pho2_eta_output,                    "pho2_eta/F");
   outTree[ii] -> Branch("pho2_phi",                &pho2_phi_output,                    "pho2_phi/F");
   outTree[ii] -> Branch("pho2_mass",               &pho2_mass_output,                   "pho2_mass/F");
   outTree[ii] -> Branch("pho2_r9",                 &pho2_r9_output,                     "pho2_r9/F");
   outTree[ii] -> Branch("pho2_sieie",              &pho2_sieie_output,                  "pho2_sieie/F");
   outTree[ii] -> Branch("pho2_hoe",                &pho2_hoe_output,                    "pho2_hoe/F");
   outTree[ii] -> Branch("jet1_pt",                 &jet1_pt_output,                     "jet1_pt/F");
   outTree[ii] -> Branch("jet1_e",                  &jet1_e_output,                      "jet1_e/F");
   outTree[ii] -> Branch("jet1_eta",                &jet1_eta_output,                    "jet1_eta/F");
   outTree[ii] -> Branch("jet1_phi",                &jet1_phi_output,                    "jet1_phi/F");
   outTree[ii] -> Branch("jet1_mass",               &jet1_mass_output,                   "jet1_mass/F");
   outTree[ii] -> Branch("jet1_csvBtag",            &jet1_csvBtag_output,                "jet1_csvBtag/F");
   outTree[ii] -> Branch("jet1_btagSF_M",           &jet1_btagSF_M_output,               "jet1_btagSF_M/F");
   outTree[ii] -> Branch("jet1_btagSFErrorUp_M",    &jet1_btagSFErrorUp_M_output,        "jet1_btagSFErrorUp_M/F");
   outTree[ii] -> Branch("jet1_btagSFErrorDown_M",  &jet1_btagSFErrorDown_M_output,      "jet1_btagSFErrorDown_M/F");
   outTree[ii] -> Branch("jet1_btagEff_M",          &jet1_btagEff_M_output,              "jet1_btagEff_M/F");
   outTree[ii] -> Branch("jet1_btagEffError_M",     &jet1_btagEffError_M_output,         "jet1_btagEffError_M/F");
   outTree[ii] -> Branch("jet1_flavour",            &jet1_flavour_output,                "jet1_flavour/I");
   outTree[ii] -> Branch("jet1_betaStarClassic",    &jet1_betaStarClassic_output,        "jet1_betaStarClassic/F");
   outTree[ii] -> Branch("jet1_dR2Mean",            &jet1_dR2Mean_output,                "jet1_dR2Mean/F");
   outTree[ii] -> Branch("jet2_pt",                 &jet2_pt_output,                     "jet2_pt/F");
   outTree[ii] -> Branch("jet2_e",                  &jet2_e_output,                      "jet2_e/F");
   outTree[ii] -> Branch("jet2_eta",                &jet2_eta_output,                    "jet2_eta/F");
   outTree[ii] -> Branch("jet2_phi",                &jet2_phi_output,                    "jet2_phi/F");
   outTree[ii] -> Branch("jet2_mass",               &jet2_mass_output,                   "jet2_mass/F");
   outTree[ii] -> Branch("jet2_csvBtag",            &jet2_csvBtag_output,                "jet2_csvBtag/F");
   outTree[ii] -> Branch("jet2_btagSF_M",           &jet2_btagSF_M_output,               "jet2_btagSF_M/F");
   outTree[ii] -> Branch("jet2_btagSFErrorUp_M",    &jet2_btagSFErrorUp_M_output,        "jet2_btagSFErrorUp_M/F");
   outTree[ii] -> Branch("jet2_btagSFErrorDown_M",  &jet2_btagSFErrorDown_M_output,      "jet2_btagSFErrorDown_M/F");
   outTree[ii] -> Branch("jet2_btagEff_M",          &jet2_btagEff_M_output,              "jet2_btagEff_M/F");
   outTree[ii] -> Branch("jet2_btagEffError_M",     &jet2_btagEffError_M_output,         "jet2_btagEffError_M/F");
   outTree[ii] -> Branch("jet2_flavour",            &jet2_flavour_output,                "jet2_flavour/I");
   outTree[ii] -> Branch("jet2_betaStarClassic",    &jet2_betaStarClassic_output,        "jet2_betaStarClassic/F");
   outTree[ii] -> Branch("jet2_dR2Mean",            &jet2_dR2Mean_output,                "jet2_dR2Mean/F");
   outTree[ii] -> Branch("jj_pt",                   &jj_pt_output,                       "jj_pt/F");
   outTree[ii] -> Branch("jj_e",                    &jj_e_output,                        "jj_e/F"); 
   outTree[ii] -> Branch("jj_phi",                  &jj_phi_output,                      "jj_phi/F");
   outTree[ii] -> Branch("jj_eta",                  &jj_eta_output,                      "jj_eta/F");
   outTree[ii] -> Branch("jj_mass",                 &jj_mass_output,                     "jj_mass/F");
   outTree[ii] -> Branch("jj_DR",                   &jj_DR_output,                       "jj_DR/F");
   outTree[ii] -> Branch("gg_pt",                   &gg_pt_output,                       "gg_pt/F");
   outTree[ii] -> Branch("gg_e",                    &gg_e_output,                        "gg_e/F");
   outTree[ii] -> Branch("gg_phi",                  &gg_phi_output,                      "gg_phi/F");
   outTree[ii] -> Branch("gg_eta",                  &gg_eta_output,                      "gg_eta/F");
   outTree[ii] -> Branch("gg_mass",                 &gg_mass_output,                     "gg_mass/F");
   outTree[ii] -> Branch("gg_DR",                   &gg_DR_output,                       "gg_DR/F");
   outTree[ii] -> Branch("costhetastar",            &costhetastar_output,                "costhetastar/F");
   outTree[ii] -> Branch("minDRgj",                 &minDRgj_output,                     "minDRgj/F");
   outTree[ii] -> Branch("ggjj_pt",                 &ggjj_pt_output,                     "ggjj_pt/F");
   outTree[ii] -> Branch("ggjj_e",                  &ggjj_e_output,                      "ggjj_e/F");
   outTree[ii] -> Branch("ggjj_phi",                &ggjj_phi_output,                    "ggjj_phi/F");
   outTree[ii] -> Branch("ggjj_eta",                &ggjj_eta_output,                    "ggjj_eta/F");
   outTree[ii] -> Branch("ggjj_mass",               &ggjj_mass_output,                   "ggjj_mass/F");
   outTree[ii] -> Branch("ggjj_DR",                 &ggjj_DR_output,                     "ggjj_DR/F");
  }
   
  TLorentzVector Pho1P4;
  TLorentzVector Pho2P4;
  TLorentzVector DiPhoP4;
  TLorentzVector jet1P4;
  TLorentzVector jet2P4;
  TLorentzVector DijetP4;
  TLorentzVector DiPhoDijetP4;
  
  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%100000==0) std::cout<<"--- Reading file:" << ntu[nn]->GetName() << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);
           
          event_output = event;
          category_output = category;
          pho1_pt_output = pho1_pt;
          pho1_e_output = pho1_e;
          pho1_phi_output = pho1_phi;
          pho1_eta_output = pho1_eta;
          pho1_mass_output = pho1_mass;
          pho1_sieie_output = pho1_sieie;
          pho1_hoe_output = pho1_hoe;
          pho1_r9_output = pho1_r9;
          pho2_pt_output = pho2_pt;
          pho2_e_output = pho2_e;
          pho2_phi_output = pho2_phi;
          pho2_eta_output = pho2_eta;
          pho2_mass_output = pho2_mass;
          pho2_r9_output = pho2_r9;
          pho2_sieie_output = pho2_sieie;
          pho2_hoe_output = pho2_hoe;
          jet1_pt_output = jet1_pt;
          jet1_e_output = jet1_e;
          jet1_phi_output = jet1_phi;
          jet1_eta_output = jet1_eta;
          jet1_mass_output = jet1_mass;
          jet1_csvBtag_output = jet1_csvBtag;
          jet1_btagSF_M_output = jet1_btagSF_M;
          jet1_btagSFErrorUp_M_output = jet1_btagSFErrorUp_M;
          jet1_btagSFErrorDown_M_output = jet1_btagSFErrorDown_M;
          jet1_btagEff_M_output = jet1_btagEff_M;
          jet1_btagEffError_M_output = jet1_btagEffError_M;
          jet1_flavour_output = jet1_flavour;
          jet1_betaStarClassic_output = jet1_betaStarClassic;
          jet2_pt_output = jet2_pt;
          jet2_e_output = jet2_e;
          jet2_phi_output = jet2_phi;
          jet2_eta_output = jet2_eta;
          jet2_mass_output = jet2_mass;
          jet2_csvBtag_output = jet2_csvBtag;
          jet2_btagSF_M_output = jet2_btagSF_M;
          jet2_btagSFErrorUp_M_output = jet2_btagSFErrorUp_M;
          jet2_btagSFErrorDown_M_output = jet2_btagSFErrorDown_M;
          jet2_btagEff_M_output = jet2_btagEff_M;
          jet2_btagEffError_M_output = jet2_btagEffError_M;
          jet2_flavour_output = jet2_flavour;
          jet2_betaStarClassic_output = jet2_betaStarClassic;

          Pho1P4.SetPtEtaPhiE(pho1_pt,pho1_eta,pho1_phi,pho1_e);
          Pho2P4.SetPtEtaPhiE(pho2_pt,pho2_eta,pho2_phi,pho2_e);
          DiPhoP4 = Pho1P4 + Pho2P4;
          jet1P4.SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2P4.SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          DijetP4 = jet1P4 + jet2P4;
          DiPhoDijetP4 = DijetP4 + DiPhoP4;
          TLorentzVector Hjj_Rstar(DijetP4);
          Hjj_Rstar.Boost(-DiPhoDijetP4.BoostVector());

          jj_pt_output = jj_pt;
          jj_e_output = jj_e;
          jj_phi_output = jj_phi;
          jj_eta_output = jj_eta;
          jj_mass_output = jj_mass;
          jj_DR_output = jj_DR;
          gg_pt_output = gg_pt;
          gg_e_output = gg_e;
          gg_phi_output = gg_phi;
          gg_eta_output = gg_eta;
          gg_mass_output = gg_mass;
          gg_DR_output = Pho1P4.DeltaR(Pho2P4);
          costhetastar_output = costhetastar;
          minDRgj_output = minDRgj;
          ggjj_pt_output = ggjj_pt;
          ggjj_phi_output = ggjj_phi;
          ggjj_eta_output = ggjj_eta;
          ggjj_mass_output = ggjj_mass;
          ggjj_DR_output = DijetP4.DeltaR(DiPhoP4);

          evweight_output = evweight;
          float weightBtagSF_output = eventWeight_2jets("medium", jet1_btagSF_M, jet2_btagSF_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          if(fabs(weightBtagSF_output) < 100.) evweight_output = weightBtagSF_output*evweight_output;

          outTree[0]->Fill();
         
      }
  }


  for(int ientry = 0; ientry < ntu_ref->GetEntries(); ientry++){
      if(ientry%100000==0) std::cout<<"--- Reading Ref_file entry = "<< ientry <<std::endl;
          ntu_ref->GetEntry(ientry);
           
          evweight = evweight*0.264/100.; //BR HbbHgg;
          
          event_output = event;
          category_output = category;
          pho1_pt_output = pho1_pt;
          pho1_e_output = pho1_e;
          pho1_phi_output = pho1_phi;
          pho1_eta_output = pho1_eta;
          pho1_mass_output = pho1_mass;
          pho1_sieie_output = pho1_sieie;
          pho1_hoe_output = pho1_hoe;
          pho1_r9_output = pho1_r9;
          pho2_pt_output = pho2_pt;
          pho2_e_output = pho2_e;
          pho2_phi_output = pho2_phi;
          pho2_eta_output = pho2_eta;
          pho2_mass_output = pho2_mass;
          pho2_r9_output = pho2_r9;
          pho2_sieie_output = pho2_sieie;
          pho2_hoe_output = pho2_hoe;
          jet1_pt_output = jet1_pt;
          jet1_e_output = jet1_e;
          jet1_phi_output = jet1_phi;
          jet1_eta_output = jet1_eta;
          jet1_mass_output = jet1_mass;
          jet1_csvBtag_output = jet1_csvBtag;
          jet1_btagSF_M_output = jet1_btagSF_M;
          jet1_btagSFErrorUp_M_output = jet1_btagSFErrorUp_M;
          jet1_btagSFErrorDown_M_output = jet1_btagSFErrorDown_M;
          jet1_btagEff_M_output = jet1_btagEff_M;
          jet1_btagEffError_M_output = jet1_btagEffError_M;
          jet1_flavour_output = jet1_flavour;
          jet1_betaStarClassic_output = jet1_betaStarClassic;
          jet2_pt_output = jet2_pt;
          jet2_e_output = jet2_e;
          jet2_phi_output = jet2_phi;
          jet2_eta_output = jet2_eta;
          jet2_mass_output = jet2_mass;
          jet2_csvBtag_output = jet2_csvBtag;
          jet2_btagSF_M_output = jet2_btagSF_M;
          jet2_btagSFErrorUp_M_output = jet2_btagSFErrorUp_M;
          jet2_btagSFErrorDown_M_output = jet2_btagSFErrorDown_M;
          jet2_btagEff_M_output = jet2_btagEff_M;
          jet2_btagEffError_M_output = jet2_btagEffError_M;
          jet2_flavour_output = jet2_flavour;
          jet2_betaStarClassic_output = jet2_betaStarClassic;

          Pho1P4.SetPtEtaPhiE(pho1_pt,pho1_eta,pho1_phi,pho1_e);
          Pho2P4.SetPtEtaPhiE(pho2_pt,pho2_eta,pho2_phi,pho2_e);
          DiPhoP4 = Pho1P4 + Pho2P4;
          jet1P4.SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2P4.SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          DijetP4 = jet1P4 + jet2P4;
          DiPhoDijetP4 = DijetP4 + DiPhoP4;
          TLorentzVector Hjj_Rstar(DijetP4);
          Hjj_Rstar.Boost(-DiPhoDijetP4.BoostVector());

          jj_pt_output = jj_pt;
          jj_e_output = jj_e;
          jj_phi_output = jj_phi;
          jj_eta_output = jj_eta;
          jj_mass_output = jj_mass;
          jj_DR_output = jj_DR;
          gg_pt_output = gg_pt;
          gg_e_output = gg_e;
          gg_phi_output = gg_phi;
          gg_eta_output = gg_eta;
          gg_mass_output = gg_mass;
          gg_DR_output = Pho1P4.DeltaR(Pho2P4);
          costhetastar_output = costhetastar;
          minDRgj_output = minDRgj;
          ggjj_pt_output = ggjj_pt;
          ggjj_phi_output = ggjj_phi;
          ggjj_eta_output = ggjj_eta;
          ggjj_mass_output = ggjj_mass;
          ggjj_DR_output = DijetP4.DeltaR(DiPhoP4);

          evweight_output = evweight;
          float weightBtagSF_output = eventWeight_2jets("medium", jet1_btagSF_M, jet2_btagSF_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          if(fabs(weightBtagSF_output) < 100.) evweight_output = weightBtagSF_output*evweight_output;

          outTree[1]->Fill();
  }

  for(int ientry = 0; ientry < ntu_data->GetEntries(); ientry++){
      if(ientry%100000==0) std::cout<<"--- Reading data_file entry = "<< ientry <<std::endl;
          ntu_data->GetEntry(ientry);
           
          event_output = event;
          category_output = category;
          pho1_pt_output = pho1_pt;
          pho1_e_output = pho1_e;
          pho1_phi_output = pho1_phi;
          pho1_eta_output = pho1_eta;
          pho1_mass_output = pho1_mass;
          pho1_sieie_output = pho1_sieie;
          pho1_hoe_output = pho1_hoe;
          pho1_r9_output = pho1_r9;
          pho2_pt_output = pho2_pt;
          pho2_e_output = pho2_e;
          pho2_phi_output = pho2_phi;
          pho2_eta_output = pho2_eta;
          pho2_mass_output = pho2_mass;
          pho2_r9_output = pho2_r9;
          pho2_sieie_output = pho2_sieie;
          pho2_hoe_output = pho2_hoe;
          jet1_pt_output = jet1_pt;
          jet1_e_output = jet1_e;
          jet1_phi_output = jet1_phi;
          jet1_eta_output = jet1_eta;
          jet1_mass_output = jet1_mass;
          jet1_csvBtag_output = jet1_csvBtag;
          jet1_btagSF_M_output = jet1_btagSF_M;
          jet1_btagSFErrorUp_M_output = jet1_btagSFErrorUp_M;
          jet1_btagSFErrorDown_M_output = jet1_btagSFErrorDown_M;
          jet1_btagEff_M_output = jet1_btagEff_M;
          jet1_btagEffError_M_output = jet1_btagEffError_M;
          jet1_flavour_output = jet1_flavour;
          jet1_betaStarClassic_output = jet1_betaStarClassic;
          jet2_pt_output = jet2_pt;
          jet2_e_output = jet2_e;
          jet2_phi_output = jet2_phi;
          jet2_eta_output = jet2_eta;
          jet2_mass_output = jet2_mass;
          jet2_csvBtag_output = jet2_csvBtag;
          jet2_btagSF_M_output = jet2_btagSF_M;
          jet2_btagSFErrorUp_M_output = jet2_btagSFErrorUp_M;
          jet2_btagSFErrorDown_M_output = jet2_btagSFErrorDown_M;
          jet2_btagEff_M_output = jet2_btagEff_M;
          jet2_btagEffError_M_output = jet2_btagEffError_M;
          jet2_flavour_output = jet2_flavour;
          jet2_betaStarClassic_output = jet2_betaStarClassic;

          Pho1P4.SetPtEtaPhiE(pho1_pt,pho1_eta,pho1_phi,pho1_e);
          Pho2P4.SetPtEtaPhiE(pho2_pt,pho2_eta,pho2_phi,pho2_e);
          DiPhoP4 = Pho1P4 + Pho2P4;
          jet1P4.SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2P4.SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          DijetP4 = jet1P4 + jet2P4;
          DiPhoDijetP4 = DijetP4 + DiPhoP4;
          TLorentzVector Hjj_Rstar(DijetP4);
          Hjj_Rstar.Boost(-DiPhoDijetP4.BoostVector());

          jj_pt_output = jj_pt;
          jj_e_output = jj_e;
          jj_phi_output = jj_phi;
          jj_eta_output = jj_eta;
          jj_mass_output = jj_mass;
          jj_DR_output = jj_DR;
          gg_pt_output = gg_pt;
          gg_e_output = gg_e;
          gg_phi_output = gg_phi;
          gg_eta_output = gg_eta;
          gg_mass_output = gg_mass;
          gg_DR_output = Pho1P4.DeltaR(Pho2P4);
          costhetastar_output = costhetastar;
          minDRgj_output = minDRgj;
          ggjj_pt_output = ggjj_pt;
          ggjj_phi_output = ggjj_phi;
          ggjj_eta_output = ggjj_eta;
          ggjj_mass_output = ggjj_mass;
          ggjj_DR_output = DijetP4.DeltaR(DiPhoP4);

          evweight_output = evweight;
          float weightBtagSF_output = eventWeight_2jets("medium", jet1_btagSF_M, jet2_btagSF_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          if(fabs(weightBtagSF_output) < 100.) evweight_output = weightBtagSF_output*evweight_output;

          outTree[2]->Fill();

  }

  TFile* outFile;
  outFile = new TFile("output/TMVA_Tree.root","RECREATE");
  outFile->cd();
  outTree[0]->Write();
  outTree[1]->Write();
  outTree[2]->Write();
  outFile->Close();
  
}



