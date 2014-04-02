
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"

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

void compareHistos(TH1F* histo_DA, TH1F* histo_MC, std::string xTitle, std::string label_DA, std::string label_MC , std::string outputDir, std::string doCUT, double& integral, double& integral_Ref);
void draw2DHisto(TH2F* h1, std::string xTitle, std::string yTitle, std::string label, std::string outputDir, std::string doCUT);
void drawNumEvents(TH1D* h1,std::string yTitle,bool& isLog, std::string label,std::string Name,std::string outputDir,std::map<int,std::string> cutString);
void compareHistosSB(TH1F* histo_DA, TH1F* histo_MC, TH1F* histo_MCbkg,std::string xTitle, std::string label_DA, std::string label_MC,std::string label_MCbkg , std::string outputDir, std::string doCUT, double& integral_data, double& integral_Ref, double& integral, std::string CAT);

int main(int argc, char** argv)
{
  
  // Set style options
  setTDRStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  //gStyle->SetOptStat(0000); 
  gStyle->SetOptFit(1); 
   
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
  std::string inputCutList = gConfigParser -> readStringOption("Input::inputCutList");

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

  char cutName[500];
  FILE *f_cuts;
  f_cuts = fopen(inputCutList.c_str(),"r");
  int pos_cut = 0;
  std::map<int,std::string> cutString;
  std::map<int,std::string> CutMap;
  std::map<int,std::vector<std::string> > CutVec;
  cutString[0] = "";
  
  while(fscanf(f_cuts,"%s \n", cutName) !=EOF ){

    std::string CUTNAME = std::string(cutName);
    if(CUTNAME.find("#") != std::string::npos) continue;
    std::cout << "\nReading input cuts: " << cutName << std::endl;

    pos_cut++;

    cutString[pos_cut] = cutString[pos_cut-1] +"_"+cutName;
    CutMap[pos_cut] = std::string(cutName);
  }  
  
  int max = 1;
  for(unsigned int ii = 1; ii < cutString.size(); ii++){
      for(int jj = 1; jj <= max; jj++){
          CutVec[ii].push_back(CutMap[jj]);
      }
      max++;
  }

  for(unsigned int ii = 1; ii < cutString.size(); ii++)
      cutString[ii] = "_cuts" + cutString[ii];

  std::map<int,std::string> catString;
  catString[0] = "";
  catString[1] = "_1btag";
  catString[2] = "_2btag";

  int           category;
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
  float         jet1_betaStarClassic;
  float         jet1_dR2Mean;
  float         jet2_pt;
  float         jet2_e;
  float         jet2_phi;
  float         jet2_eta;
  float         jet2_mass;
  float         jet2_csvBtag;
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
   TBranch        *b_jet1_betaStarClassic;   //!
   TBranch        *b_jet1_dR2Mean;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet2_e;   //!
   TBranch        *b_jet2_phi;   //!
   TBranch        *b_jet2_eta;   //!
   TBranch        *b_jet2_mass;   //!
   TBranch        *b_jet2_csvBtag;   //!
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
   ntu[ii]->SetBranchAddress("jet1_betaStarClassic", &jet1_betaStarClassic, &b_jet1_betaStarClassic);
   ntu[ii]->SetBranchAddress("jet1_dR2Mean", &jet1_dR2Mean, &b_jet1_dR2Mean);
   ntu[ii]->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   ntu[ii]->SetBranchAddress("jet2_e", &jet2_e, &b_jet2_e);
   ntu[ii]->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   ntu[ii]->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   ntu[ii]->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   ntu[ii]->SetBranchAddress("jet2_csvBtag", &jet2_csvBtag, &b_jet2_csvBtag);
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
   ntu_ref->SetBranchAddress("jet1_betaStarClassic", &jet1_betaStarClassic, &b_jet1_betaStarClassic);
   ntu_ref->SetBranchAddress("jet1_dR2Mean", &jet1_dR2Mean, &b_jet1_dR2Mean);
   ntu_ref->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   ntu_ref->SetBranchAddress("jet2_e", &jet2_e, &b_jet2_e);
   ntu_ref->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   ntu_ref->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   ntu_ref->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   ntu_ref->SetBranchAddress("jet2_csvBtag", &jet2_csvBtag, &b_jet2_csvBtag);
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
   ntu_data->SetBranchAddress("jet1_betaStarClassic", &jet1_betaStarClassic, &b_jet1_betaStarClassic);
   ntu_data->SetBranchAddress("jet1_dR2Mean", &jet1_dR2Mean, &b_jet1_dR2Mean);
   ntu_data->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   ntu_data->SetBranchAddress("jet2_e", &jet2_e, &b_jet2_e);
   ntu_data->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   ntu_data->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   ntu_data->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   ntu_data->SetBranchAddress("jet2_csvBtag", &jet2_csvBtag, &b_jet2_csvBtag);
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

  TH1D* h_numberEvents = new TH1D("h_numberEvents","h_numberEvents",cutString.size(),0,cutString.size());
  TH1D* h_ref_numberEvents = new TH1D("h_ref_numberEvents","h_ref_numberEvents",cutString.size(),0,cutString.size());
  TH1D* h_FigOfMerit = new TH1D("h_FigOfMerit","h_FigOfMerit",cutString.size(),0,cutString.size());

  TH1D* h_numberEvents_1btag = new TH1D("h_numberEvents_1btag","h_numberEvents_1btag",cutString.size(),0,cutString.size());
  TH1D* h_ref_numberEvents_1btag = new TH1D("h_ref_numberEvents_1btag","h_ref_numberEvents_1btag",cutString.size(),0,cutString.size());
  TH1D* h_FigOfMerit_1btag = new TH1D("h_FigOfMerit_1btag","h_FigOfMerit_1btag",cutString.size(),0,cutString.size());

  TH1D* h_numberEvents_2btag = new TH1D("h_numberEvents_2btag","h_numberEvents_2btag",cutString.size(),0,cutString.size());
  TH1D* h_ref_numberEvents_2btag = new TH1D("h_ref_numberEvents_2btag","h_ref_numberEvents_2btag",cutString.size(),0,cutString.size());
  TH1D* h_FigOfMerit_2btag = new TH1D("h_FigOfMerit_2btag","h_FigOfMerit_2btag",cutString.size(),0,cutString.size());
   
  std::map<int,TH1F*> h_NEvents;
  std::map<int,TH1F*> h_NEvents_1btag;
  std::map<int,TH1F*> h_NEvents_2btag;
  std::map<int,TH1F*> h_NEvents_SB;
  std::map<int,TH1F*> h_NEvents_1btag_SB;
  std::map<int,TH1F*> h_NEvents_2btag_SB;

  std::map<int,std::map<int,TH1F*> > h_PhoPt; 
  std::map<int,std::map<int,TH1F*> > h_PhoEta;
  std::map<int,std::map<int,TH1F*> > h_Pho1Eta;
  std::map<int,std::map<int,TH1F*> > h_Pho2Eta;
  std::map<int,std::map<int,TH1F*> > h_PhoPhi; 
  std::map<int,std::map<int,TH1F*> > h_PhoR9; 
  std::map<int,std::map<int,TH1F*> > h_PhoHoE; 
  std::map<int,std::map<int,TH1F*> > h_PhoSieie; 
  std::map<int,std::map<int,TH1F*> > h_DiPhoPt;
  std::map<int,TH1F*> h_DiPhoPt_SB;
  std::map<int,TH1F*> h_DiPhoPt_1btag_SB;
  std::map<int,TH1F*> h_DiPhoPt_2btag_SB;
  std::map<int,std::map<int,TH1F*> > h_DiPhoEta;
  std::map<int,std::map<int,TH1F*> > h_DiPhoPhi;
  std::map<int,std::map<int,TH1F*> > h_DiPhoDeltaR;
  std::map<int,std::map<int,TH1F*> > h_DiPhoInvMass;
  std::map<int,std::map<int,TH1F*> > h_DiPhoCosthetastar;
  std::map<int,std::map<int,TH1F*> > h_jetPt;
  std::map<int,std::map<int,TH1F*> > h_jetEta;
  std::map<int,std::map<int,TH1F*> > h_jetPhi;
  std::map<int,std::map<int,TH1F*> > h_jetcsvBtag;
  std::map<int,std::map<int,TH1F*> > h_jetbetaStarClassic;
  std::map<int,std::map<int,TH1F*> > h_jetdR2Mean;
  std::map<int,std::map<int,TH1F*> > h_DijetPt;
  std::map<int,std::map<int,TH1F*> > h_DijetEta;
  std::map<int,std::map<int,TH1F*> > h_DijetPhi;
  std::map<int,std::map<int,TH1F*> > h_DijetDeltaR;
  std::map<int,std::map<int,TH1F*> > h_DijetInvMass;
  std::map<int,std::map<int,TH1F*> > h_DijetCosthetastar;
  std::map<int,std::map<int,TH1F*> > h_minPhojetDeltaR;
  std::map<int,std::map<int,TH1F*> > h_DiPhoDijetDeltaR;
  std::map<int,std::map<int,TH1F*> > h_DiPhoDijetPt;
  std::map<int,std::map<int,TH1F*> > h_DiPhoDijetEta;
  std::map<int,std::map<int,TH1F*> > h_DiPhoDijetPhi;
  std::map<int,std::map<int,TH1F*> > h_DiPhoDijetInvMass;
  std::map<int,std::map<int,TH2F*> > h2_Pho1Pt_vs_jet1Pt;
  std::map<int,std::map<int,TH2F*> > h2_Pho2Pt_vs_jet2Pt;
  std::map<int,std::map<int,TH2F*> > h2_Pho1Eta_vs_jet1Eta;
  std::map<int,std::map<int,TH2F*> > h2_Pho2Eta_vs_jet2Eta;
  std::map<int,std::map<int,TH2F*> > h2_Pho1Phi_vs_jet1Phi;
  std::map<int,std::map<int,TH2F*> > h2_Pho2Phi_vs_jet2Phi;
  std::map<int,std::map<int,TH2F*> > h2_DiPhoInvMass_vs_DijetInvMass;
  std::map<int,std::map<int,TH2F*> > h2_DiPhoPt_vs_DijetPt;
  std::map<int,std::map<int,TH2F*> > h2_DiPhoEta_vs_DijetEta;
  std::map<int,std::map<int,TH2F*> > h2_DiPhoPhi_vs_DijetPhi;
  std::map<int,std::map<int,TH2F*> > h2_DiPhoDeltaR_vs_DijetDeltaR;
  std::map<int,std::map<int,TH2F*> > h2_DiPhoCosthetastar_vs_DijetCosthetastar;


  for(unsigned int ii = 0; ii < cutString.size(); ii++){
   h_NEvents[ii] = new TH1F((std::string("h_NEvents")+cutString[ii]).c_str(),(std::string("h_NEvents")+cutString[ii]).c_str(),100000000,0,100000000);
   h_NEvents_1btag[ii] = new TH1F((std::string("h_NEvents_1btag")+cutString[ii]).c_str(),(std::string("h_NEvents_1btag")+cutString[ii]).c_str(),100000000,0,100000000);
   h_NEvents_2btag[ii] = new TH1F((std::string("h_NEvents_2btag")+cutString[ii]).c_str(),(std::string("h_NEvents_2btag")+cutString[ii]).c_str(),100000000,0,100000000);
   h_NEvents_SB[ii] = new TH1F((std::string("h_NEvents_SB")+cutString[ii]).c_str(),(std::string("h_NEvents_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_NEvents_1btag_SB[ii] = new TH1F((std::string("h_NEvents_1btag_SB")+cutString[ii]).c_str(),(std::string("h_NEvents_1btag")+cutString[ii]).c_str(),100000000,0,100000000);
   h_NEvents_2btag_SB[ii] = new TH1F((std::string("h_NEvents_2btag_SB")+cutString[ii]).c_str(),(std::string("h_NEvents_2btag_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_DiPhoPt_SB[ii] = new TH1F((std::string("h_DiPhoPt_SB")+cutString[ii]).c_str(),(std::string("h_DiPhoPt_SB")+cutString[ii]).c_str(),120,0.,600.);  
   h_DiPhoPt_1btag_SB[ii] = new TH1F((std::string("h_DiPhoPt_1btag_SB")+cutString[ii]).c_str(),(std::string("h_DiPhoPt_1btag_SB")+cutString[ii]).c_str(),120,0.,600.);
   h_DiPhoPt_2btag_SB[ii] = new TH1F((std::string("h_DiPhoPt_2btag_SB")+cutString[ii]).c_str(),(std::string("h_DiPhoPt_2btag_SB")+cutString[ii]).c_str(),120,0.,600.);
  }
 
 for(unsigned int ii = 0; ii < 3; ii++){
  for(unsigned int jj = 0; jj < cutString.size(); jj++){
   h_PhoPt[ii][jj] = new TH1F((std::string("h_PhoPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_PhoPt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.);
   h_PhoEta[ii][jj] = new TH1F((std::string("h_PhoEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_PhoEta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_Pho1Eta[ii][jj] = new TH1F((std::string("h_Pho1Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_Pho1Eta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_Pho2Eta[ii][jj] = new TH1F((std::string("h_Pho2Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_Pho2Eta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_PhoPhi[ii][jj] = new TH1F((std::string("h_PhoPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_PhoPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_PhoR9[ii][jj] = new TH1F((std::string("h_PhoR9")+cutString[jj]+catString[ii]).c_str(),(std::string("h_PhoR9")+cutString[jj]+catString[ii]).c_str(),50,0.,1.);
   h_PhoHoE[ii][jj] = new TH1F((std::string("h_PhoHoE")+cutString[jj]+catString[ii]).c_str(),(std::string("h_PhoHoE")+cutString[jj]+catString[ii]).c_str(),50,0.,0.05);
   h_PhoSieie[ii][jj] = new TH1F((std::string("h_PhoSieie")+cutString[jj]+catString[ii]).c_str(),(std::string("h_PhoSieie")+cutString[jj]+catString[ii]).c_str(),50,0.,0.07);
   h_DiPhoPt[ii][jj] = new TH1F((std::string("h_DiPhoPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.);
   h_DiPhoEta[ii][jj] = new TH1F((std::string("h_DiPhoEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoEta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_DiPhoPhi[ii][jj] = new TH1F((std::string("h_DiPhoPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_DiPhoDeltaR[ii][jj] = new TH1F((std::string("h_DiPhoDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoDeltaR")+cutString[jj]+catString[ii]).c_str(),60,0.,6.);
   h_DiPhoInvMass[ii][jj] = new TH1F((std::string("h_DiPhoInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoInvMass")+cutString[jj]+catString[ii]).c_str(),50,100.,180.);
   h_DiPhoCosthetastar[ii][jj] = new TH1F((std::string("h_DiPhoCosthetastar")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoCosthetastar")+cutString[jj]+catString[ii]).c_str(),60,-1.5,1.5);
   h_jetPt[ii][jj] = new TH1F((std::string("h_jetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_jetPt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.);
   h_jetEta[ii][jj] = new TH1F((std::string("h_jetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_jetEta")+cutString[jj]+catString[ii]).c_str(),80,-4,4);
   h_jetPhi[ii][jj] = new TH1F((std::string("h_jetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_jetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_jetcsvBtag[ii][jj] = new TH1F((std::string("h_jetcsvBtag")+cutString[jj]+catString[ii]).c_str(),(std::string("h_jetcsvBtag")+cutString[jj]+catString[ii]).c_str(),100,0.,1.);
   h_jetbetaStarClassic[ii][jj] = new TH1F((std::string("h_jetbetaStarClassic")+cutString[jj]+catString[ii]).c_str(),(std::string("h_jetbetaStarClassic")+cutString[jj]+catString[ii]).c_str(),100,0.,1.);
   h_jetdR2Mean[ii][jj] = new TH1F((std::string("h_jetdR2Mean")+cutString[jj]+catString[ii]).c_str(),(std::string("h_jetdR2Mean")+cutString[jj]+catString[ii]).c_str(),100,0.,0.1);
   h_DijetPt[ii][jj] = new TH1F((std::string("h_DijetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DijetPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.);
   h_DijetEta[ii][jj] = new TH1F((std::string("h_DijetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DijetEta")+cutString[jj]+catString[ii]).c_str(),80,-4,4);
   h_DijetPhi[ii][jj] = new TH1F((std::string("h_DijetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DijetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_DijetDeltaR[ii][jj] = new TH1F((std::string("h_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),70,0.,7.);
   h_DijetInvMass[ii][jj] = new TH1F((std::string("h_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),60,0.,300.);
   h_DijetCosthetastar[ii][jj] = new TH1F((std::string("h_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),60,-1.5,1.5);
   h_minPhojetDeltaR[ii][jj] = new TH1F((std::string("h_minPhojetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_minPhojetDeltaR")+cutString[jj]+catString[ii]).c_str(),80,0.,4.);
   h_DiPhoDijetDeltaR[ii][jj] = new TH1F((std::string("h_DiPhoDijetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoDijetDeltaR")+cutString[jj]+catString[ii]).c_str(),90,0.,9.);
   h_DiPhoDijetPt[ii][jj] = new TH1F((std::string("h_DiPhoDijetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoDijetPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.);
   h_DiPhoDijetEta[ii][jj] = new TH1F((std::string("h_DiPhoDijetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoDijetEta")+cutString[jj]+catString[ii]).c_str(),100,-5,5);
   h_DiPhoDijetPhi[ii][jj] = new TH1F((std::string("h_DiPhoDijetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoDijetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_DiPhoDijetInvMass[ii][jj] = new TH1F((std::string("h_DiPhoDijetInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h_DiPhoDijetInvMass")+cutString[jj]+catString[ii]).c_str(),200,0.,1200.); 
   h2_Pho1Pt_vs_jet1Pt[ii][jj] = new TH2F((std::string("h2_Pho1Pt_vs_jet1Pt")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_Pho1Pt_vs_jet1Pt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.,100,0.,500.);
   h2_Pho2Pt_vs_jet2Pt[ii][jj] = new TH2F((std::string("h2_Pho2Pt_vs_jet2Pt")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_Pho2Pt_vs_jet2Pt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.,100,0.,500.);
   h2_Pho1Eta_vs_jet1Eta[ii][jj] = new TH2F((std::string("h2_Pho1Eta_vs_jet1Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_Pho1Eta_vs_jet1Eta")+cutString[jj]+catString[ii]).c_str(),80,-4,4,80,-4,4);
   h2_Pho2Eta_vs_jet2Eta[ii][jj] = new TH2F((std::string("h2_Pho2Eta_vs_jet2Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_Pho2Eta_vs_jet2Eta")+cutString[jj]+catString[ii]).c_str(),80,-4,4,80,-4,4);
   h2_Pho1Phi_vs_jet1Phi[ii][jj] = new TH2F((std::string("h2_Pho1Phi_vs_jet1Phi")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_Pho1Phi_vs_jet1Phi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_Pho2Phi_vs_jet2Phi[ii][jj] = new TH2F((std::string("h2_Pho2Phi_vs_jet2Phi")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_Pho2Phi_vs_jet2Phi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_DiPhoInvMass_vs_DijetInvMass[ii][jj] = new TH2F((std::string("h2_DiPhoInvMass_vs_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_DiPhoInvMass_vs_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),50,100.,180.,60,0.,300.);
   h2_DiPhoPt_vs_DijetPt[ii][jj] = new TH2F((std::string("h2_DiPhoPt_vs_DijetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_DiPhoPt_vs_DijetPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.,120,0.,600.);
   h2_DiPhoEta_vs_DijetEta[ii][jj] = new TH2F((std::string("h2_DiPhoEta_vs_DijetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_DiPhoEta_vs_DijetEta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_DiPhoPhi_vs_DijetPhi[ii][jj] = new TH2F((std::string("h2_DiPhoPhi_vs_DijetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_DiPhoPhi_vs_DijetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_DiPhoDeltaR_vs_DijetDeltaR[ii][jj] = new TH2F((std::string("h2_DiPhoDeltaR_vs_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_DiPhoDeltaR_vs_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),60,0.,6.,60,0.,6.);
   h2_DiPhoCosthetastar_vs_DijetCosthetastar[ii][jj] = new TH2F((std::string("h2_DiPhoCosthetastar_vs_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_DiPhoCosthetastar_vs_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),60,-1.5,1.5,60,-1.5,1.5);
  }
 }
  
  std::map<int,TH1F*> h_ref_NEvents;
  std::map<int,TH1F*> h_ref_NEvents_1btag;
  std::map<int,TH1F*> h_ref_NEvents_2btag;
  std::map<int,TH1F*> h_ref_NEvents_SB;
  std::map<int,TH1F*> h_ref_NEvents_1btag_SB;
  std::map<int,TH1F*> h_ref_NEvents_2btag_SB;

  std::map<int,std::map<int,TH1F*> > h_ref_PhoPt; 
  std::map<int,std::map<int,TH1F*> > h_ref_PhoEta;
  std::map<int,std::map<int,TH1F*> > h_ref_Pho1Eta;
  std::map<int,std::map<int,TH1F*> > h_ref_Pho2Eta;
  std::map<int,std::map<int,TH1F*> > h_ref_PhoPhi; 
  std::map<int,std::map<int,TH1F*> > h_ref_PhoR9; 
  std::map<int,std::map<int,TH1F*> > h_ref_PhoHoE; 
  std::map<int,std::map<int,TH1F*> > h_ref_PhoSieie; 
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoPt;
  std::map<int,TH1F*> h_ref_DiPhoPt_SB;
  std::map<int,TH1F*> h_ref_DiPhoPt_1btag_SB;
  std::map<int,TH1F*> h_ref_DiPhoPt_2btag_SB;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoEta;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoPhi;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoDeltaR;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoInvMass;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoCosthetastar;
  std::map<int,std::map<int,TH1F*> > h_ref_jetPt;
  std::map<int,std::map<int,TH1F*> > h_ref_jetEta;
  std::map<int,std::map<int,TH1F*> > h_ref_jetPhi;
  std::map<int,std::map<int,TH1F*> > h_ref_jetcsvBtag;
  std::map<int,std::map<int,TH1F*> > h_ref_jetbetaStarClassic;
  std::map<int,std::map<int,TH1F*> > h_ref_jetdR2Mean;
  std::map<int,std::map<int,TH1F*> > h_ref_DijetPt;
  std::map<int,std::map<int,TH1F*> > h_ref_DijetEta;
  std::map<int,std::map<int,TH1F*> > h_ref_DijetPhi;
  std::map<int,std::map<int,TH1F*> > h_ref_DijetDeltaR;
  std::map<int,std::map<int,TH1F*> > h_ref_DijetInvMass;
  std::map<int,std::map<int,TH1F*> > h_ref_DijetCosthetastar;
  std::map<int,std::map<int,TH1F*> > h_ref_minPhojetDeltaR;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoDijetDeltaR;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoDijetPt;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoDijetEta;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoDijetPhi;
  std::map<int,std::map<int,TH1F*> > h_ref_DiPhoDijetInvMass;
  std::map<int,std::map<int,TH2F*> > h2_ref_Pho1Pt_vs_jet1Pt;
  std::map<int,std::map<int,TH2F*> > h2_ref_Pho2Pt_vs_jet2Pt;
  std::map<int,std::map<int,TH2F*> > h2_ref_Pho1Eta_vs_jet1Eta;
  std::map<int,std::map<int,TH2F*> > h2_ref_Pho2Eta_vs_jet2Eta;
  std::map<int,std::map<int,TH2F*> > h2_ref_Pho1Phi_vs_jet1Phi;
  std::map<int,std::map<int,TH2F*> > h2_ref_Pho2Phi_vs_jet2Phi;
  std::map<int,std::map<int,TH2F*> > h2_ref_DiPhoInvMass_vs_DijetInvMass;
  std::map<int,std::map<int,TH2F*> > h2_ref_DiPhoPt_vs_DijetPt;
  std::map<int,std::map<int,TH2F*> > h2_ref_DiPhoEta_vs_DijetEta;
  std::map<int,std::map<int,TH2F*> > h2_ref_DiPhoPhi_vs_DijetPhi;
  std::map<int,std::map<int,TH2F*> > h2_ref_DiPhoDeltaR_vs_DijetDeltaR;
  std::map<int,std::map<int,TH2F*> > h2_ref_DiPhoCosthetastar_vs_DijetCosthetastar;

  for(unsigned int ii = 0; ii < cutString.size(); ii++){
   h_ref_NEvents[ii] = new TH1F((std::string("h_ref_NEvents")+cutString[ii]).c_str(),(std::string("h_ref_NEvents")+cutString[ii]).c_str(),100000000,0,100000000);
   h_ref_NEvents_1btag[ii] = new TH1F((std::string("h_ref_NEvents_1btag")+cutString[ii]).c_str(),(std::string("h_ref_NEvents_1btag")+cutString[ii]).c_str(),100000000,0,100000000);
   h_ref_NEvents_2btag[ii] = new TH1F((std::string("h_ref_NEvents_2btag")+cutString[ii]).c_str(),(std::string("h_ref_NEvents_2btag")+cutString[ii]).c_str(),100000000,0,100000000);
   h_ref_NEvents_SB[ii] = new TH1F((std::string("h_ref_NEvents_SB")+cutString[ii]).c_str(),(std::string("h_ref_NEvents_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_ref_NEvents_1btag_SB[ii] = new TH1F((std::string("h_ref_NEvents_1btag_SB")+cutString[ii]).c_str(),(std::string("h_ref_NEvents_1btag")+cutString[ii]).c_str(),100000000,0,100000000);
   h_ref_NEvents_2btag_SB[ii] = new TH1F((std::string("h_ref_NEvents_2btag_SB")+cutString[ii]).c_str(),(std::string("h_ref_NEvents_2btag_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_ref_DiPhoPt_SB[ii] = new TH1F((std::string("h_ref_DiPhoPt_SB")+cutString[ii]).c_str(),(std::string("h_ref_DiPhoPt_SB")+cutString[ii]).c_str(),120,0.,600.);  
   h_ref_DiPhoPt_1btag_SB[ii] = new TH1F((std::string("h_ref_DiPhoPt_1btag_SB")+cutString[ii]).c_str(),(std::string("h_ref_DiPhoPt_1btag_SB")+cutString[ii]).c_str(),120,0.,600.);
   h_ref_DiPhoPt_2btag_SB[ii] = new TH1F((std::string("h_ref_DiPhoPt_2btag_SB")+cutString[ii]).c_str(),(std::string("h_ref_DiPhoPt_2btag_SB")+cutString[ii]).c_str(),120,0.,600.);
  }
 
 for(unsigned int ii = 0; ii < 3; ii++){
  for(unsigned int jj = 0; jj < cutString.size(); jj++){
   h_ref_PhoPt[ii][jj] = new TH1F((std::string("h_ref_PhoPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_PhoPt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.);
   h_ref_PhoEta[ii][jj] = new TH1F((std::string("h_ref_PhoEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_PhoEta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_Pho1Eta[ii][jj] = new TH1F((std::string("h_ref_Pho1Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_Pho1Eta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_Pho2Eta[ii][jj] = new TH1F((std::string("h_ref_Pho2Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_Pho2Eta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_PhoPhi[ii][jj] = new TH1F((std::string("h_ref_PhoPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_PhoPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_PhoR9[ii][jj] = new TH1F((std::string("h_ref_PhoR9")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_PhoR9")+cutString[jj]+catString[ii]).c_str(),50,0.,1.);
   h_ref_PhoHoE[ii][jj] = new TH1F((std::string("h_ref_PhoHoE")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_PhoHoE")+cutString[jj]+catString[ii]).c_str(),50,0.,0.05);
   h_ref_PhoSieie[ii][jj] = new TH1F((std::string("h_ref_PhoSieie")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_PhoSieie")+cutString[jj]+catString[ii]).c_str(),50,0.,0.07);
   h_ref_DiPhoPt[ii][jj] = new TH1F((std::string("h_ref_DiPhoPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.);
   h_ref_DiPhoEta[ii][jj] = new TH1F((std::string("h_ref_DiPhoEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoEta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_DiPhoPhi[ii][jj] = new TH1F((std::string("h_ref_DiPhoPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_DiPhoDeltaR[ii][jj] = new TH1F((std::string("h_ref_DiPhoDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoDeltaR")+cutString[jj]+catString[ii]).c_str(),60,0.,6.);
   h_ref_DiPhoInvMass[ii][jj] = new TH1F((std::string("h_ref_DiPhoInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoInvMass")+cutString[jj]+catString[ii]).c_str(),50,100.,180.);
   h_ref_DiPhoCosthetastar[ii][jj] = new TH1F((std::string("h_ref_DiPhoCosthetastar")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoCosthetastar")+cutString[jj]+catString[ii]).c_str(),60,-1.5,1.5);
   h_ref_jetPt[ii][jj] = new TH1F((std::string("h_ref_jetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_jetPt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.);
   h_ref_jetEta[ii][jj] = new TH1F((std::string("h_ref_jetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_jetEta")+cutString[jj]+catString[ii]).c_str(),80,-4,4);
   h_ref_jetPhi[ii][jj] = new TH1F((std::string("h_ref_jetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_jetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_jetcsvBtag[ii][jj] = new TH1F((std::string("h_ref_jetcsvBtag")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_jetcsvBtag")+cutString[jj]+catString[ii]).c_str(),100,0.,1.);
   h_ref_jetbetaStarClassic[ii][jj] = new TH1F((std::string("h_ref_jetbetaStarClassic")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_jetbetaStarClassic")+cutString[jj]+catString[ii]).c_str(),100,0.,1.);
   h_ref_jetdR2Mean[ii][jj] = new TH1F((std::string("h_ref_jetdR2Mean")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_jetdR2Mean")+cutString[jj]+catString[ii]).c_str(),100,0.,0.1);
   h_ref_DijetPt[ii][jj] = new TH1F((std::string("h_ref_DijetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DijetPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.);
   h_ref_DijetEta[ii][jj] = new TH1F((std::string("h_ref_DijetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DijetEta")+cutString[jj]+catString[ii]).c_str(),80,-4,4);
   h_ref_DijetPhi[ii][jj] = new TH1F((std::string("h_ref_DijetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DijetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_DijetDeltaR[ii][jj] = new TH1F((std::string("h_ref_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),70,0.,7.);
   h_ref_DijetInvMass[ii][jj] = new TH1F((std::string("h_ref_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),60,0.,300.);
   h_ref_DijetCosthetastar[ii][jj] = new TH1F((std::string("h_ref_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),60,-1.5,1.5);
   h_ref_minPhojetDeltaR[ii][jj] = new TH1F((std::string("h_ref_minPhojetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_minPhojetDeltaR")+cutString[jj]+catString[ii]).c_str(),80,0.,4.);
   h_ref_DiPhoDijetDeltaR[ii][jj] = new TH1F((std::string("h_ref_DiPhoDijetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoDijetDeltaR")+cutString[jj]+catString[ii]).c_str(),90,0.,9.);
   h_ref_DiPhoDijetPt[ii][jj] = new TH1F((std::string("h_ref_DiPhoDijetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoDijetPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.);
   h_ref_DiPhoDijetEta[ii][jj] = new TH1F((std::string("h_ref_DiPhoDijetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoDijetEta")+cutString[jj]+catString[ii]).c_str(),100,-5,5);
   h_ref_DiPhoDijetPhi[ii][jj] = new TH1F((std::string("h_ref_DiPhoDijetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoDijetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5);
   h_ref_DiPhoDijetInvMass[ii][jj] = new TH1F((std::string("h_ref_DiPhoDijetInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h_ref_DiPhoDijetInvMass")+cutString[jj]+catString[ii]).c_str(),200,0.,1200.); 
   h2_ref_Pho1Pt_vs_jet1Pt[ii][jj] = new TH2F((std::string("h2_ref_Pho1Pt_vs_jet1Pt")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_Pho1Pt_vs_jet1Pt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.,100,0.,500.);
   h2_ref_Pho2Pt_vs_jet2Pt[ii][jj] = new TH2F((std::string("h2_ref_Pho2Pt_vs_jet2Pt")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_Pho2Pt_vs_jet2Pt")+cutString[jj]+catString[ii]).c_str(),100,0.,500.,100,0.,500.);
   h2_ref_Pho1Eta_vs_jet1Eta[ii][jj] = new TH2F((std::string("h2_ref_Pho1Eta_vs_jet1Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_Pho1Eta_vs_jet1Eta")+cutString[jj]+catString[ii]).c_str(),80,-4,4,80,-4,4);
   h2_ref_Pho2Eta_vs_jet2Eta[ii][jj] = new TH2F((std::string("h2_ref_Pho2Eta_vs_jet2Eta")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_Pho2Eta_vs_jet2Eta")+cutString[jj]+catString[ii]).c_str(),80,-4,4,80,-4,4);
   h2_ref_Pho1Phi_vs_jet1Phi[ii][jj] = new TH2F((std::string("h2_ref_Pho1Phi_vs_jet1Phi")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_Pho1Phi_vs_jet1Phi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_ref_Pho2Phi_vs_jet2Phi[ii][jj] = new TH2F((std::string("h2_ref_Pho2Phi_vs_jet2Phi")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_Pho2Phi_vs_jet2Phi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_ref_DiPhoInvMass_vs_DijetInvMass[ii][jj] = new TH2F((std::string("h2_ref_DiPhoInvMass_vs_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_DiPhoInvMass_vs_DijetInvMass")+cutString[jj]+catString[ii]).c_str(),50,100.,180.,60,0.,300.);
   h2_ref_DiPhoPt_vs_DijetPt[ii][jj] = new TH2F((std::string("h2_ref_DiPhoPt_vs_DijetPt")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_DiPhoPt_vs_DijetPt")+cutString[jj]+catString[ii]).c_str(),120,0.,600.,120,0.,600.);
   h2_ref_DiPhoEta_vs_DijetEta[ii][jj] = new TH2F((std::string("h2_ref_DiPhoEta_vs_DijetEta")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_DiPhoEta_vs_DijetEta")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_ref_DiPhoPhi_vs_DijetPhi[ii][jj] = new TH2F((std::string("h2_ref_DiPhoPhi_vs_DijetPhi")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_DiPhoPhi_vs_DijetPhi")+cutString[jj]+catString[ii]).c_str(),70,-3.5,3.5,70,-3.5,3.5);
   h2_ref_DiPhoDeltaR_vs_DijetDeltaR[ii][jj] = new TH2F((std::string("h2_ref_DiPhoDeltaR_vs_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_DiPhoDeltaR_vs_DijetDeltaR")+cutString[jj]+catString[ii]).c_str(),60,0.,6.,60,0.,6.);
   h2_ref_DiPhoCosthetastar_vs_DijetCosthetastar[ii][jj] = new TH2F((std::string("h2_ref_DiPhoCosthetastar_vs_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),(std::string("h2_ref_DiPhoCosthetastar_vs_DijetCosthetastar")+cutString[jj]+catString[ii]).c_str(),60,-1.5,1.5,60,-1.5,1.5);
  }
 }
  
  std::map<int,TH1F*> h_data_NEvents_SB;
  std::map<int,TH1F*> h_data_NEvents_1btag_SB;
  std::map<int,TH1F*> h_data_NEvents_2btag_SB;
  std::map<int,TH1F*> h_data_DiPhoPt_SB;
  std::map<int,TH1F*> h_data_DiPhoPt_1btag_SB;
  std::map<int,TH1F*> h_data_DiPhoPt_2btag_SB;

  for(unsigned int ii = 0; ii < cutString.size(); ii++){
   h_data_NEvents_SB[ii] = new TH1F((std::string("h_data_NEvents_SB")+cutString[ii]).c_str(),(std::string("h_data_NEvents_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_data_NEvents_1btag_SB[ii] = new TH1F((std::string("h_data_NEvents_1btag_SB")+cutString[ii]).c_str(),(std::string("h_data_NEvents_1btag_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_data_NEvents_2btag_SB[ii] = new TH1F((std::string("h_data_NEvents_2btag_SB")+cutString[ii]).c_str(),(std::string("h_data_NEvents_2btag_SB")+cutString[ii]).c_str(),100000000,0,100000000);
   h_data_DiPhoPt_SB[ii] = new TH1F((std::string("h_data_DiPhoPt_SB")+cutString[ii]).c_str(),(std::string("h_data_DiPhoPt_SB")+cutString[ii]).c_str(),120,0.,600.);
   h_data_DiPhoPt_1btag_SB[ii] = new TH1F((std::string("h_data_DiPhoPt_1btag_SB")+cutString[ii]).c_str(),(std::string("h_data_DiPhoPt_1btag_SB")+cutString[ii]).c_str(),120,0.,600.);
   h_data_DiPhoPt_2btag_SB[ii] = new TH1F((std::string("h_data_DiPhoPt_2btag_SB")+cutString[ii]).c_str(),(std::string("h_data_DiPhoPt_2btag_SB")+cutString[ii]).c_str(),120,0.,600.);
  }

  TLorentzVector Pho1P4;
  TLorentzVector Pho2P4;
  TLorentzVector DiPhoP4;
  TLorentzVector jet1P4;
  TLorentzVector jet2P4;
  TLorentzVector DijetP4;
  TLorentzVector DiPhoDijetP4;

  std::map<std::string,bool> typeCut;

  int ievent[cutString.size()];
  int ievent_1btag[cutString.size()];
  int ievent_2btag[cutString.size()];
  int ievent_SB[cutString.size()];
  int ievent_1btag_SB[cutString.size()];
  int ievent_2btag_SB[cutString.size()];
  for(unsigned int ii = 0; ii<cutString.size();ii++){
      ievent[ii] = 0;
      ievent_1btag[ii] = 0;
      ievent_2btag[ii] = 0;
      ievent_SB[ii] = 0;
      ievent_1btag_SB[ii] = 0;
      ievent_2btag_SB[ii] = 0;
  }

  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%100000==0) std::cout<<"--- Reading file:" << ntu[nn]->GetName() << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

          Pho1P4.SetPtEtaPhiE(pho1_pt,pho1_eta,pho1_phi,pho1_e);
          Pho2P4.SetPtEtaPhiE(pho2_pt,pho2_eta,pho2_phi,pho2_e);
          DiPhoP4 = Pho1P4 + Pho2P4;
          jet1P4.SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2P4.SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          DijetP4 = jet1P4 + jet2P4;
          DiPhoDijetP4 = DijetP4 + DiPhoP4;
          TLorentzVector Hjj_Rstar(DijetP4);
          Hjj_Rstar.Boost(-DiPhoDijetP4.BoostVector());

          if(cutString[cutString.size()-1].find("ggMass") != std::string::npos && (DiPhoP4.M() >= 120. && DiPhoP4.M() <= 130.)) typeCut[std::string("ggMass")] = true;
          else typeCut[std::string("ggMass")] = false;

          if(cutString[cutString.size()-1].find("jjMass") != std::string::npos && (DijetP4.M() >= 85. && DijetP4.M() <= 155.)) typeCut[std::string("jjMass")] = true;
          else typeCut[std::string("jjMass")] = false;

          if(cutString[cutString.size()-1].find("ggCosthetastar") != std::string::npos && fabs(costhetastar) <= 0.65) typeCut[std::string("ggCosthetastar")] = true;
          else typeCut[std::string("ggCosthetastar")] = false;

          if(cutString[cutString.size()-1].find("jjCosthetastar") != std::string::npos && fabs(costhetastar) <= 0.65) typeCut[std::string("jjCosthetastar")] = true;
          else typeCut[std::string("jjCosthetastar")] = false;

          if(cutString[cutString.size()-1].find("ggDR") != std::string::npos && Pho1P4.DeltaR(Pho2P4) < 2.) typeCut[std::string("ggDR")] = true;
          else typeCut[std::string("ggDR")] = false;

          if(cutString[cutString.size()-1].find("jjDR") != std::string::npos && jj_DR < 2.) typeCut[std::string("jjDR")] = true;
          else typeCut[std::string("jjDR")] = false;

          if(cutString[cutString.size()-1].find("gEta") != std::string::npos && (fabs(pho1_eta) < 1.5 && fabs(pho2_eta) < 1.5)) typeCut[std::string("gEta")] = true;
          else typeCut[std::string("gEta")] = false;
          
          for(unsigned int ii = 0; ii < cutString.size(); ii++){
             
              bool isGood = true;
              if(ii != 0)
                 for(unsigned int jj = 0; jj < CutVec[ii].size(); jj++){
                     
                     if(typeCut[CutVec[ii].at(jj)] == false) isGood = false;
                 }

              if(isGood == false) continue;

              if(DiPhoP4.M() >= 120. && DiPhoP4.M() <= 130.){

                 ievent[ii]++;
                 h_NEvents[ii]->Fill(ievent[ii],evweight);

                 if(category == 1) ievent_1btag[ii]++;
                 if(category == 1) h_NEvents_1btag[ii]->Fill(ievent_1btag[ii],evweight);

                 if(category == 2) ievent_2btag[ii]++;
                 if(category == 2) h_NEvents_2btag[ii]->Fill(ievent_2btag[ii],evweight);
              }

              if(DiPhoP4.M() < 120. || DiPhoP4.M() > 130.){

                 ievent_SB[ii]++;
                 h_NEvents_SB[ii]->Fill(ievent_SB[ii],evweight);
                 h_DiPhoPt_SB[ii]->Fill(gg_pt,evweight);

                 if(category == 1) ievent_1btag_SB[ii]++;
                 if(category == 1) h_NEvents_1btag_SB[ii]->Fill(ievent_1btag_SB[ii],evweight);
                 if(category == 1) h_DiPhoPt_1btag_SB[ii]->Fill(gg_pt,evweight);

                 if(category == 2) ievent_2btag_SB[ii]++;
                 if(category == 2) h_NEvents_2btag_SB[ii]->Fill(ievent_2btag_SB[ii],evweight);
                 if(category == 2) h_DiPhoPt_2btag_SB[ii]->Fill(gg_pt,evweight);
              }
              
              for(int jj = 0; jj < 3; jj++){
              
               if(jj == 1 && category != 1) continue;
               if(jj == 2 && category != 2) continue;
               
               h_PhoPt[jj][ii]->Fill(pho1_pt,evweight);       h_PhoPt[jj][ii]->Fill(pho2_pt,evweight);
               h_PhoEta[jj][ii]->Fill(pho1_eta,evweight);     h_PhoEta[jj][ii]->Fill(pho2_eta,evweight);
               h_Pho1Eta[jj][ii]->Fill(pho1_eta,evweight);    h_Pho2Eta[jj][ii]->Fill(pho2_eta,evweight);
               h_PhoPhi[jj][ii]->Fill(pho1_phi,evweight);     h_PhoPhi[jj][ii]->Fill(pho2_phi,evweight);
               h_PhoR9[jj][ii]->Fill(pho1_r9,evweight);       h_PhoR9[jj][ii]->Fill(pho2_r9,evweight);
               h_PhoHoE[jj][ii]->Fill(pho1_hoe,evweight);     h_PhoHoE[jj][ii]->Fill(pho2_hoe,evweight);
               h_PhoSieie[jj][ii]->Fill(pho1_sieie,evweight); h_PhoSieie[jj][ii]->Fill(pho2_sieie,evweight);

               h_DiPhoPt[jj][ii]->Fill(gg_pt,evweight);
               h_DiPhoEta[jj][ii]->Fill(gg_eta,evweight);
               h_DiPhoPhi[jj][ii]->Fill(gg_phi,evweight);
               h_DiPhoDeltaR[jj][ii]->Fill(Pho1P4.DeltaR(Pho2P4),evweight);
               h_DiPhoInvMass[jj][ii]->Fill(gg_mass,evweight);
               h_DiPhoCosthetastar[jj][ii]->Fill(costhetastar,evweight);
      
               h_jetPt[jj][ii]->Fill(jet1_pt,evweight);                           h_jetPt[jj][ii]->Fill(jet2_pt,evweight);
               h_jetEta[jj][ii]->Fill(jet1_eta,evweight);                         h_jetEta[jj][ii]->Fill(jet2_eta,evweight);
               h_jetPhi[jj][ii]->Fill(jet1_phi,evweight);                         h_jetPhi[jj][ii]->Fill(jet2_phi,evweight);
               h_jetcsvBtag[jj][ii]->Fill(jet1_csvBtag,evweight);                 h_jetcsvBtag[jj][ii]->Fill(jet2_csvBtag,evweight);
               h_jetbetaStarClassic[jj][ii]->Fill(jet1_betaStarClassic,evweight); h_jetbetaStarClassic[jj][ii]->Fill(jet2_betaStarClassic,evweight);
               h_jetdR2Mean[jj][ii]->Fill(jet1_dR2Mean,evweight);                 h_jetdR2Mean[jj][ii]->Fill(jet2_dR2Mean,evweight);
          
               h_DijetPt[jj][ii]->Fill(jj_pt,evweight);
               h_DijetEta[jj][ii]->Fill(jj_eta,evweight);
               h_DijetPhi[jj][ii]->Fill(jj_phi,evweight);
               h_DijetDeltaR[jj][ii]->Fill(jj_DR,evweight);
               h_DijetInvMass[jj][ii]->Fill(DijetP4.M(),evweight);
               h_DijetCosthetastar[jj][ii]->Fill(Hjj_Rstar.CosTheta(),evweight);

               h_minPhojetDeltaR[jj][ii]->Fill(minDRgj,evweight); 

               h_DiPhoDijetDeltaR[jj][ii]->Fill(DijetP4.DeltaR(DiPhoP4),evweight);
               h_DiPhoDijetPt[jj][ii]->Fill(DiPhoDijetP4.Pt(),evweight);
               h_DiPhoDijetEta[jj][ii]->Fill(DiPhoDijetP4.Eta(),evweight);
               h_DiPhoDijetPhi[jj][ii]->Fill(DiPhoDijetP4.Phi(),evweight);
               h_DiPhoDijetInvMass[jj][ii]->Fill(DiPhoDijetP4.M(),evweight);

               if(evweight > 5) continue;

               h2_Pho1Pt_vs_jet1Pt[jj][ii]->Fill(pho1_pt,jet1_pt,evweight);
               h2_Pho2Pt_vs_jet2Pt[jj][ii]->Fill(pho2_pt,jet2_pt,evweight);
               h2_Pho1Eta_vs_jet1Eta[jj][ii]->Fill(pho1_eta,jet1_eta,evweight);
               h2_Pho2Eta_vs_jet2Eta[jj][ii]->Fill(pho2_eta,jet2_eta,evweight);
               h2_Pho1Phi_vs_jet1Phi[jj][ii]->Fill(pho1_phi,jet1_phi,evweight);
               h2_Pho2Phi_vs_jet2Phi[jj][ii]->Fill(pho2_phi,jet2_phi,evweight);
               h2_DiPhoInvMass_vs_DijetInvMass[jj][ii]->Fill(gg_mass,DijetP4.M(),evweight);
               h2_DiPhoPt_vs_DijetPt[jj][ii]->Fill(gg_pt,jj_pt,evweight);
               h2_DiPhoEta_vs_DijetEta[jj][ii]->Fill(gg_eta,jj_eta,evweight);
               h2_DiPhoPhi_vs_DijetPhi[jj][ii]->Fill(gg_phi,jj_phi,evweight);
               h2_DiPhoDeltaR_vs_DijetDeltaR[jj][ii]->Fill(Pho1P4.DeltaR(Pho2P4),jj_DR,evweight);
               h2_DiPhoCosthetastar_vs_DijetCosthetastar[jj][ii]->Fill(costhetastar,Hjj_Rstar.CosTheta(),evweight);
              }
          }
      }
  }

  for(unsigned int ii = 0; ii<cutString.size();ii++){
      ievent[ii] = 0;
      ievent_1btag[ii] = 0;
      ievent_2btag[ii] = 0;
      ievent_SB[ii] = 0;
      ievent_1btag_SB[ii] = 0;
      ievent_2btag_SB[ii] = 0;
  }
 
  for(int ientry = 0; ientry < ntu_ref->GetEntries(); ientry++){
      if(ientry%100000==0) std::cout<<"--- Reading Ref_file entry = "<< ientry <<std::endl;
          ntu_ref->GetEntry(ientry);
           
          evweight = evweight*0.264/100.; //BR HbbHgg;
          
          Pho1P4.SetPtEtaPhiE(pho1_pt,pho1_eta,pho1_phi,pho1_e);
          Pho2P4.SetPtEtaPhiE(pho2_pt,pho2_eta,pho2_phi,pho2_e);
          DiPhoP4 = Pho1P4 + Pho2P4;
          jet1P4.SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2P4.SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          DijetP4 = jet1P4 + jet2P4;
          DiPhoDijetP4 = DijetP4 + DiPhoP4;
          TLorentzVector Hjj_Rstar(DijetP4);
          Hjj_Rstar.Boost(-DiPhoDijetP4.BoostVector());

          if(cutString[cutString.size()-1].find("ggMass") != std::string::npos && (DiPhoP4.M() >= 120. && DiPhoP4.M() <= 130.)) typeCut[std::string("ggMass")] = true;
          else typeCut[std::string("ggMass")] = false;

          if(cutString[cutString.size()-1].find("jjMass") != std::string::npos && (DijetP4.M() >= 85. && DijetP4.M() <= 155.)) typeCut[std::string("jjMass")] = true;
          else typeCut[std::string("jjMass")] = false;

          if(cutString[cutString.size()-1].find("ggCosthetastar") != std::string::npos && fabs(costhetastar) <= 0.65) typeCut[std::string("ggCosthetastar")] = true;
          else typeCut[std::string("ggCosthetastar")] = false;

          if(cutString[cutString.size()-1].find("jjCosthetastar") != std::string::npos && fabs(costhetastar) <= 0.65) typeCut[std::string("jjCosthetastar")] = true;
          else typeCut[std::string("jjCosthetastar")] = false;

          if(cutString[cutString.size()-1].find("ggDR") != std::string::npos && Pho1P4.DeltaR(Pho2P4) < 2.) typeCut[std::string("ggDR")] = true;
          else typeCut[std::string("ggDR")] = false;

          if(cutString[cutString.size()-1].find("jjDR") != std::string::npos && jj_DR < 2.) typeCut[std::string("jjDR")] = true;
          else typeCut[std::string("jjDR")] = false;

          if(cutString[cutString.size()-1].find("gEta") != std::string::npos && (fabs(pho1_eta) < 1.5 && fabs(pho2_eta) < 1.5)) typeCut[std::string("gEta")] = true;
          else typeCut[std::string("gEta")] = false;
          
          for(unsigned int ii = 0; ii < cutString.size(); ii++){
             
              bool isGood = true;
              if(ii != 0)
                 for(unsigned int jj = 0; jj < CutVec[ii].size(); jj++){
                     
                     if(typeCut[CutVec[ii].at(jj)] == false) isGood = false;
                 }

              if(isGood == false) continue;

              if(DiPhoP4.M() >= 120. && DiPhoP4.M() <= 130.){

                 ievent[ii]++;
                 h_ref_NEvents[ii]->Fill(ievent[ii],evweight);

                 if(category == 1) ievent_1btag[ii]++;
                 if(category == 1) h_ref_NEvents_1btag[ii]->Fill(ievent_1btag[ii],evweight);

                 if(category == 2) ievent_2btag[ii]++;
                 if(category == 2) h_ref_NEvents_2btag[ii]->Fill(ievent_2btag[ii],evweight);
              }

              if(DiPhoP4.M() < 120. || DiPhoP4.M() > 130.){

                 ievent_SB[ii]++;
                 h_ref_NEvents_SB[ii]->Fill(ievent_SB[ii],evweight);
                 h_ref_DiPhoPt_SB[ii]->Fill(gg_pt,evweight);

                 if(category == 1) ievent_1btag_SB[ii]++;
                 if(category == 1) h_ref_NEvents_1btag_SB[ii]->Fill(ievent_1btag_SB[ii],evweight);
                 if(category == 1) h_ref_DiPhoPt_1btag_SB[ii]->Fill(gg_pt,evweight);

                 if(category == 2) ievent_2btag_SB[ii]++;
                 if(category == 2) h_ref_NEvents_2btag_SB[ii]->Fill(ievent_2btag_SB[ii],evweight);
                 if(category == 2) h_ref_DiPhoPt_2btag_SB[ii]->Fill(gg_pt,evweight);
              }

              for(int jj = 0; jj < 3; jj++){
              
               if(jj == 1 && category != 1) continue;
               if(jj == 2 && category != 2) continue;
               
               h_ref_PhoPt[jj][ii]->Fill(pho1_pt,evweight);       h_ref_PhoPt[jj][ii]->Fill(pho2_pt,evweight);
               h_ref_PhoEta[jj][ii]->Fill(pho1_eta,evweight);     h_ref_PhoEta[jj][ii]->Fill(pho2_eta,evweight);
               h_ref_Pho1Eta[jj][ii]->Fill(pho1_eta,evweight);    h_ref_Pho2Eta[jj][ii]->Fill(pho2_eta,evweight);
               h_ref_PhoPhi[jj][ii]->Fill(pho1_phi,evweight);     h_ref_PhoPhi[jj][ii]->Fill(pho2_phi,evweight);
               h_ref_PhoR9[jj][ii]->Fill(pho1_r9,evweight);       h_ref_PhoR9[jj][ii]->Fill(pho2_r9,evweight);
               h_ref_PhoHoE[jj][ii]->Fill(pho1_hoe,evweight);     h_ref_PhoHoE[jj][ii]->Fill(pho2_hoe,evweight);
               h_ref_PhoSieie[jj][ii]->Fill(pho1_sieie,evweight); h_ref_PhoSieie[jj][ii]->Fill(pho2_sieie,evweight);

               h_ref_DiPhoPt[jj][ii]->Fill(gg_pt,evweight);
               h_ref_DiPhoEta[jj][ii]->Fill(gg_eta,evweight);
               h_ref_DiPhoPhi[jj][ii]->Fill(gg_phi,evweight);
               h_ref_DiPhoDeltaR[jj][ii]->Fill(Pho1P4.DeltaR(Pho2P4),evweight);
               h_ref_DiPhoInvMass[jj][ii]->Fill(gg_mass,evweight);
               h_ref_DiPhoCosthetastar[jj][ii]->Fill(costhetastar,evweight);
      
               h_ref_jetPt[jj][ii]->Fill(jet1_pt,evweight);                           h_ref_jetPt[jj][ii]->Fill(jet2_pt,evweight);
               h_ref_jetEta[jj][ii]->Fill(jet1_eta,evweight);                         h_ref_jetEta[jj][ii]->Fill(jet2_eta,evweight);
               h_ref_jetPhi[jj][ii]->Fill(jet1_phi,evweight);                         h_ref_jetPhi[jj][ii]->Fill(jet2_phi,evweight);
               h_ref_jetcsvBtag[jj][ii]->Fill(jet1_csvBtag,evweight);                 h_ref_jetcsvBtag[jj][ii]->Fill(jet2_csvBtag,evweight);
               h_ref_jetbetaStarClassic[jj][ii]->Fill(jet1_betaStarClassic,evweight); h_ref_jetbetaStarClassic[jj][ii]->Fill(jet2_betaStarClassic,evweight);
               h_ref_jetdR2Mean[jj][ii]->Fill(jet1_dR2Mean,evweight);                 h_ref_jetdR2Mean[jj][ii]->Fill(jet2_dR2Mean,evweight);
          
               h_ref_DijetPt[jj][ii]->Fill(jj_pt,evweight);
               h_ref_DijetEta[jj][ii]->Fill(jj_eta,evweight);
               h_ref_DijetPhi[jj][ii]->Fill(jj_phi,evweight);
               h_ref_DijetDeltaR[jj][ii]->Fill(jj_DR,evweight);
               h_ref_DijetInvMass[jj][ii]->Fill(DijetP4.M(),evweight);
               h_ref_DijetCosthetastar[jj][ii]->Fill(Hjj_Rstar.CosTheta(),evweight);

               h_ref_minPhojetDeltaR[jj][ii]->Fill(minDRgj,evweight); 

               h_ref_DiPhoDijetDeltaR[jj][ii]->Fill(DijetP4.DeltaR(DiPhoP4),evweight);
               h_ref_DiPhoDijetPt[jj][ii]->Fill(DiPhoDijetP4.Pt(),evweight);
               h_ref_DiPhoDijetEta[jj][ii]->Fill(DiPhoDijetP4.Eta(),evweight);
               h_ref_DiPhoDijetPhi[jj][ii]->Fill(DiPhoDijetP4.Phi(),evweight);
               h_ref_DiPhoDijetInvMass[jj][ii]->Fill(DiPhoDijetP4.M(),evweight);

               if(evweight > 5) continue;

               h2_ref_Pho1Pt_vs_jet1Pt[jj][ii]->Fill(pho1_pt,jet1_pt,evweight);
               h2_ref_Pho2Pt_vs_jet2Pt[jj][ii]->Fill(pho2_pt,jet2_pt,evweight);
               h2_ref_Pho1Eta_vs_jet1Eta[jj][ii]->Fill(pho1_eta,jet1_eta,evweight);
               h2_ref_Pho2Eta_vs_jet2Eta[jj][ii]->Fill(pho2_eta,jet2_eta,evweight);
               h2_ref_Pho1Phi_vs_jet1Phi[jj][ii]->Fill(pho1_phi,jet1_phi,evweight);
               h2_ref_Pho2Phi_vs_jet2Phi[jj][ii]->Fill(pho2_phi,jet2_phi,evweight);
               h2_ref_DiPhoInvMass_vs_DijetInvMass[jj][ii]->Fill(gg_mass,DijetP4.M(),evweight);
               h2_ref_DiPhoPt_vs_DijetPt[jj][ii]->Fill(gg_pt,jj_pt,evweight);
               h2_ref_DiPhoEta_vs_DijetEta[jj][ii]->Fill(gg_eta,jj_eta,evweight);
               h2_ref_DiPhoPhi_vs_DijetPhi[jj][ii]->Fill(gg_phi,jj_phi,evweight);
               h2_ref_DiPhoDeltaR_vs_DijetDeltaR[jj][ii]->Fill(Pho1P4.DeltaR(Pho2P4),jj_DR,evweight);
               h2_ref_DiPhoCosthetastar_vs_DijetCosthetastar[jj][ii]->Fill(costhetastar,Hjj_Rstar.CosTheta(),evweight);
              }
          }    

  }

  for(unsigned int ii = 0; ii<cutString.size();ii++){
      ievent[ii] = 0;
      ievent_1btag[ii] = 0;
      ievent_2btag[ii] = 0;
      ievent_SB[ii] = 0;
      ievent_1btag_SB[ii] = 0;
      ievent_2btag_SB[ii] = 0;
  }

  for(int ientry = 0; ientry < ntu_data->GetEntries(); ientry++){
      if(ientry%100000==0) std::cout<<"--- Reading data_file entry = "<< ientry <<std::endl;
          ntu_data->GetEntry(ientry);
           
          Pho1P4.SetPtEtaPhiE(pho1_pt,pho1_eta,pho1_phi,pho1_e);
          Pho2P4.SetPtEtaPhiE(pho2_pt,pho2_eta,pho2_phi,pho2_e);
          DiPhoP4 = Pho1P4 + Pho2P4;
          jet1P4.SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2P4.SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          DijetP4 = jet1P4 + jet2P4;
          DiPhoDijetP4 = DijetP4 + DiPhoP4;
          TLorentzVector Hjj_Rstar(DijetP4);
          Hjj_Rstar.Boost(-DiPhoDijetP4.BoostVector());

          if(cutString[cutString.size()-1].find("ggMass") != std::string::npos && (DiPhoP4.M() >= 120. && DiPhoP4.M() <= 130.)) typeCut[std::string("ggMass")] = true;
          else typeCut[std::string("ggMass")] = false;

          if(cutString[cutString.size()-1].find("jjMass") != std::string::npos && (DijetP4.M() >= 85. && DijetP4.M() <= 155.)) typeCut[std::string("jjMass")] = true;
          else typeCut[std::string("jjMass")] = false;

          if(cutString[cutString.size()-1].find("ggCosthetastar") != std::string::npos && fabs(costhetastar) <= 0.65) typeCut[std::string("ggCosthetastar")] = true;
          else typeCut[std::string("ggCosthetastar")] = false;

          if(cutString[cutString.size()-1].find("jjCosthetastar") != std::string::npos && fabs(costhetastar) <= 0.65) typeCut[std::string("jjCosthetastar")] = true;
          else typeCut[std::string("jjCosthetastar")] = false;

          if(cutString[cutString.size()-1].find("ggDR") != std::string::npos && Pho1P4.DeltaR(Pho2P4) < 2.) typeCut[std::string("ggDR")] = true;
          else typeCut[std::string("ggDR")] = false;

          if(cutString[cutString.size()-1].find("jjDR") != std::string::npos && jj_DR < 2.) typeCut[std::string("jjDR")] = true;
          else typeCut[std::string("jjDR")] = false;

          if(cutString[cutString.size()-1].find("gEta") != std::string::npos && (fabs(pho1_eta) < 1.5 && fabs(pho2_eta) < 1.5)) typeCut[std::string("gEta")] = true;
          else typeCut[std::string("gEta")] = false;
          
          for(unsigned int ii = 0; ii < cutString.size(); ii++){
             
              bool isGood = true;
              if(ii != 0)
                 for(unsigned int jj = 0; jj < CutVec[ii].size(); jj++){
                     
                     if(typeCut[CutVec[ii].at(jj)] == false) isGood = false;
                 }

              if(isGood == false) continue;

              if(DiPhoP4.M() >= 120. && DiPhoP4.M() <= 130.) continue;

              ievent_SB[ii]++;
              h_data_NEvents_SB[ii]->Fill(ievent_SB[ii]);
              h_data_DiPhoPt_SB[ii]->Fill(gg_pt);

              if(category == 1) ievent_1btag_SB[ii]++;
              if(category == 1) h_data_NEvents_1btag_SB[ii]->Fill(ievent_1btag_SB[ii]);
              if(category == 1) h_data_DiPhoPt_1btag_SB[ii]->Fill(gg_pt);

              if(category == 2) ievent_2btag_SB[ii]++;
              if(category == 2) h_data_NEvents_2btag_SB[ii]->Fill(ievent_2btag_SB[ii]);
              if(category == 2) h_data_DiPhoPt_2btag_SB[ii]->Fill(gg_pt); 
          }

  }
 
  for(unsigned int ii = 0; ii < cutString.size(); ii++){

      double integral = h_NEvents[ii]->Integral();
      double integral_Ref = h_ref_NEvents[ii]->Integral();
      
      double integral_1btag = h_NEvents_1btag[ii]->Integral();
      double integral_Ref_1btag = h_ref_NEvents_1btag[ii]->Integral();

      double integral_2btag = h_NEvents_2btag[ii]->Integral();
      double integral_Ref_2btag = h_ref_NEvents_2btag[ii]->Integral();

      double integral_SB = h_NEvents_SB[ii]->Integral();
      double integral_Ref_SB = h_ref_NEvents_SB[ii]->Integral();
      double integral_data_SB = h_data_NEvents_SB[ii]->Integral();
      
      double integral_1btag_SB = h_NEvents_1btag_SB[ii]->Integral();
      double integral_Ref_1btag_SB = h_ref_NEvents_1btag_SB[ii]->Integral();
      double integral_data_1btag_SB = h_data_NEvents_1btag_SB[ii]->Integral();

      double integral_2btag_SB = h_NEvents_2btag_SB[ii]->Integral();
      double integral_Ref_2btag_SB = h_ref_NEvents_2btag_SB[ii]->Integral();
      double integral_data_2btag_SB = h_data_NEvents_2btag_SB[ii]->Integral();

      h_numberEvents->SetBinContent(h_numberEvents->FindBin(ii),integral);
      h_ref_numberEvents->SetBinContent(h_ref_numberEvents->FindBin(ii),integral_Ref);
      h_FigOfMerit->SetBinContent(h_FigOfMerit->FindBin(ii),integral_Ref/sqrt(integral_Ref+integral));

      h_numberEvents_1btag->SetBinContent(h_numberEvents_1btag->FindBin(ii),integral_1btag);
      h_ref_numberEvents_1btag->SetBinContent(h_ref_numberEvents_1btag->FindBin(ii),integral_Ref_1btag);
      h_FigOfMerit_1btag->SetBinContent(h_FigOfMerit_1btag->FindBin(ii),integral_Ref_1btag/sqrt(integral_Ref_1btag+integral_1btag));
      
      h_numberEvents_2btag->SetBinContent(h_numberEvents_2btag->FindBin(ii),integral_2btag);
      h_ref_numberEvents_2btag->SetBinContent(h_ref_numberEvents_2btag->FindBin(ii),integral_Ref_2btag);
      h_FigOfMerit_2btag->SetBinContent(h_FigOfMerit_2btag->FindBin(ii),integral_Ref_2btag/sqrt(integral_Ref_2btag+integral_2btag));

      for(int jj = 0; jj < 3; jj++){
       std::string doCUT = cutString[ii]+catString[jj];
       gStyle->SetOptStat(1110); 
       //gStyle->SetOptStat(0000); 

       double INT;
       double INT_REF;

       if(jj == 0){
          INT = h_NEvents[ii]->Integral();
          INT_REF = h_ref_NEvents[ii]->Integral();
       }
       if(jj == 1){
          INT = h_NEvents_1btag[ii]->Integral();
          INT_REF = h_ref_NEvents_1btag[ii]->Integral();
       }
       if(jj == 2){
          INT = h_NEvents_2btag[ii]->Integral();
          INT_REF = h_ref_NEvents_2btag[ii]->Integral();
       }
       
       compareHistos(h_ref_PhoPt[jj][ii], h_PhoPt[jj][ii], std::string("g_pt"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_PhoEta[jj][ii], h_PhoEta[jj][ii], std::string("g_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_Pho1Eta[jj][ii], h_Pho1Eta[jj][ii], std::string("g1_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_Pho2Eta[jj][ii], h_Pho2Eta[jj][ii], std::string("g2_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_PhoPhi[jj][ii], h_PhoPhi[jj][ii], std::string("g_phi"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_PhoR9[jj][ii], h_PhoR9[jj][ii], std::string("g_r9"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_PhoHoE[jj][ii], h_PhoHoE[jj][ii], std::string("g_hoe"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_PhoSieie[jj][ii], h_PhoSieie[jj][ii], std::string("g_sieie"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);

       compareHistos(h_ref_DiPhoPt[jj][ii], h_DiPhoPt[jj][ii], std::string("gg_pt"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoEta[jj][ii], h_DiPhoEta[jj][ii], std::string("gg_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoPhi[jj][ii], h_DiPhoPhi[jj][ii], std::string("gg_phi"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoDeltaR[jj][ii], h_DiPhoDeltaR[jj][ii], std::string("gg_dR"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoInvMass[jj][ii], h_DiPhoInvMass[jj][ii], std::string("gg_Mass"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoCosthetastar[jj][ii], h_DiPhoCosthetastar[jj][ii], std::string("gg_costhetastar"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
  
       compareHistos(h_ref_jetPt[jj][ii], h_jetPt[jj][ii], std::string("j_pt"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_jetEta[jj][ii], h_jetEta[jj][ii], std::string("j_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_jetPhi[jj][ii], h_jetPhi[jj][ii], std::string("j_phi"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_jetcsvBtag[jj][ii], h_jetcsvBtag[jj][ii], std::string("j_csvBtag"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);  
       compareHistos(h_ref_jetbetaStarClassic[jj][ii], h_jetbetaStarClassic[jj][ii], std::string("j_betaStarClassic"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_jetdR2Mean[jj][ii], h_jetdR2Mean[jj][ii], std::string("j_dR2Mean"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);

       compareHistos(h_ref_DijetPt[jj][ii], h_DijetPt[jj][ii], std::string("jj_pt"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DijetEta[jj][ii], h_DijetEta[jj][ii], std::string("jj_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DijetPhi[jj][ii], h_DijetPhi[jj][ii], std::string("jj_phi"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DijetDeltaR[jj][ii], h_DijetDeltaR[jj][ii], std::string("jj_dR"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DijetInvMass[jj][ii], h_DijetInvMass[jj][ii], std::string("jj_Mass"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DijetCosthetastar[jj][ii], h_DijetCosthetastar[jj][ii], std::string("jj_costhetastar"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);

       compareHistos(h_ref_minPhojetDeltaR[jj][ii], h_minPhojetDeltaR[jj][ii], std::string("minGJ_dR"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF); 

       compareHistos(h_ref_DiPhoDijetDeltaR[jj][ii], h_DiPhoDijetDeltaR[jj][ii], std::string("ggjj_dR"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoDijetPt[jj][ii], h_DiPhoDijetPt[jj][ii], std::string("ggjj_pt"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoDijetEta[jj][ii], h_DiPhoDijetEta[jj][ii], std::string("ggjj_eta"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoDijetPhi[jj][ii], h_DiPhoDijetPhi[jj][ii], std::string("ggjj_phi"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF);
       compareHistos(h_ref_DiPhoDijetInvMass[jj][ii], h_DiPhoDijetInvMass[jj][ii], std::string("ggjj_mass"), std::string("signal"), std::string("bkg") , outputDir, doCUT, INT, INT_REF); 

       //gStyle->SetOptStat(1110); 
       gStyle->SetOptStat(0000); 
       draw2DHisto(h2_ref_Pho1Pt_vs_jet1Pt[jj][ii], std::string("g1_pt"), std::string("j1_pt"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_Pho1Pt_vs_jet1Pt[jj][ii], std::string("g1_pt"), std::string("j1_pt"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_Pho2Pt_vs_jet2Pt[jj][ii], std::string("g2_pt"), std::string("j2_pt"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_Pho2Pt_vs_jet2Pt[jj][ii], std::string("g2_pt"), std::string("j2_pt"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_Pho1Eta_vs_jet1Eta[jj][ii], std::string("g1_eta"), std::string("j1_eta"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_Pho1Eta_vs_jet1Eta[jj][ii], std::string("g1_eta"), std::string("j1_eta"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_Pho2Eta_vs_jet2Eta[jj][ii], std::string("g2_eta"), std::string("j2_eta"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_Pho2Eta_vs_jet2Eta[jj][ii], std::string("g2_eta"), std::string("j2_eta"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_Pho1Phi_vs_jet1Phi[jj][ii], std::string("g1_phi"), std::string("j1_phi"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_Pho1Phi_vs_jet1Phi[jj][ii], std::string("g1_phi"), std::string("j1_phi"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_Pho2Phi_vs_jet2Phi[jj][ii], std::string("g2_phi"), std::string("j2_phi"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_Pho2Phi_vs_jet2Phi[jj][ii], std::string("g2_phi"), std::string("j2_phi"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_DiPhoInvMass_vs_DijetInvMass[jj][ii], std::string("gg_mass"), std::string("jj_mass"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_DiPhoInvMass_vs_DijetInvMass[jj][ii], std::string("gg_mass"), std::string("jj_mass"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_DiPhoPt_vs_DijetPt[jj][ii], std::string("gg_pt"), std::string("jj_pt"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_DiPhoPt_vs_DijetPt[jj][ii], std::string("gg_pt"), std::string("jj_pt"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_DiPhoEta_vs_DijetEta[jj][ii], std::string("gg_eta"), std::string("jj_eta"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_DiPhoEta_vs_DijetEta[jj][ii], std::string("gg_eta"), std::string("jj_eta"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_DiPhoPhi_vs_DijetPhi[jj][ii], std::string("gg_phi"), std::string("jj_phi"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_DiPhoPhi_vs_DijetPhi[jj][ii], std::string("gg_phi"), std::string("jj_phi"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_DiPhoDeltaR_vs_DijetDeltaR[jj][ii], std::string("gg_dR"), std::string("jj_dR"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_DiPhoDeltaR_vs_DijetDeltaR[jj][ii], std::string("gg_dR"), std::string("jj_dR"),std::string("bkg"),outputDir,doCUT);
       draw2DHisto(h2_ref_DiPhoCosthetastar_vs_DijetCosthetastar[jj][ii], std::string("gg_costhetastar"), std::string("jj_costhetastar"),std::string("signal"),outputDir,doCUT);
       draw2DHisto(h2_DiPhoCosthetastar_vs_DijetCosthetastar[jj][ii], std::string("gg_costhetastar"), std::string("jj_costhetastar"),std::string("bkg"),outputDir,doCUT); 
      }
      
      std::string doCUT = cutString[ii];
      compareHistosSB(h_data_DiPhoPt_SB[ii],h_ref_DiPhoPt_SB[ii],h_DiPhoPt_SB[ii],std::string("gg_pt"), std::string("Data"), std::string("Signal"),std::string("Bkg") , outputDir, doCUT,integral_data_SB,integral_Ref_SB,integral_SB, std::string(""));
      compareHistosSB(h_data_DiPhoPt_1btag_SB[ii],h_ref_DiPhoPt_1btag_SB[ii],h_DiPhoPt_1btag_SB[ii],std::string("gg_pt"), std::string("Data"), std::string("Signal"),std::string("Bkg") , outputDir, doCUT,integral_data_1btag_SB,integral_Ref_1btag_SB,integral_1btag_SB,std::string("_1btag"));
      compareHistosSB(h_data_DiPhoPt_2btag_SB[ii],h_ref_DiPhoPt_2btag_SB[ii],h_DiPhoPt_2btag_SB[ii],std::string("gg_pt"), std::string("Data"), std::string("Signal"),std::string("Bkg") , outputDir, doCUT,integral_data_2btag_SB,integral_Ref_2btag_SB,integral_2btag_SB,std::string("_2btag"));
      
  }

  bool isLog = true;
  drawNumEvents(h_numberEvents,std::string("events"),isLog,std::string("bkg"),std::string("NumberEvents"),outputDir,cutString);
  drawNumEvents(h_ref_numberEvents,std::string("events"),isLog,std::string("signal"),std::string("NumberEvents"),outputDir,cutString);
  drawNumEvents(h_numberEvents_1btag,std::string("events"),isLog,std::string("bkg"),std::string("NumberEvents_1btag"),outputDir,cutString);
  drawNumEvents(h_ref_numberEvents_1btag,std::string("events"),isLog,std::string("signal"),std::string("NumberEvents_1btag"),outputDir,cutString);
  drawNumEvents(h_numberEvents_2btag,std::string("events"),isLog,std::string("bkg"),std::string("NumberEvents_2btag"),outputDir,cutString);
  drawNumEvents(h_ref_numberEvents_2btag,std::string("events"),isLog,std::string("signal"),std::string("NumberEvents_2btag"),outputDir,cutString);
  isLog = false;
  drawNumEvents(h_FigOfMerit,std::string("S/sqrt(S+B)"),isLog,std::string("bkg-signal"),std::string("FigOfMerit"),outputDir,cutString);
  drawNumEvents(h_FigOfMerit_1btag,std::string("S/sqrt(S+B)"),isLog,std::string("bkg-signal"),std::string("FigOfMerit_1btag"),outputDir,cutString);
  drawNumEvents(h_FigOfMerit_2btag,std::string("S/sqrt(S+B)"),isLog,std::string("bkg-signal"),std::string("FigOfMerit_2btag"),outputDir,cutString);

}


void compareHistos(TH1F* histo_DA, TH1F* histo_MC, std::string xTitle, std::string label_DA, std::string label_MC , std::string outputDir, std::string doCUT, double& integral, double& integral_Ref){   
    
    char intDA_buffer[50];
    sprintf(intDA_buffer,("#evt-"+label_DA+" = %1.2e").c_str(),integral_Ref);

    char intMC_buffer[50];
    sprintf(intMC_buffer,("#evt-"+label_MC+" = %1.2e").c_str(),integral);

    double FoM = integral_Ref/sqrt(integral+integral_Ref);
    char FoM_buffer[50];
    sprintf(FoM_buffer,"S/sqrt(S+B) = %f",FoM);

    histo_DA -> Scale(1./histo_DA->Integral());
    histo_MC -> Scale(1./histo_MC->Integral());

    histo_DA -> GetXaxis() -> SetLabelSize(0);
    histo_DA -> GetYaxis() -> SetLabelSize(0.04);
    histo_DA -> GetXaxis() -> SetTitleSize(0.05);
    histo_DA -> GetYaxis() -> SetTitleSize(0.05);
    histo_DA -> GetYaxis() -> SetTitleOffset(1.25);

    histo_DA -> GetXaxis() -> SetTitle(xTitle.c_str());
    histo_MC -> GetXaxis() -> SetTitle(xTitle.c_str());
    histo_DA -> GetYaxis() -> SetTitle("");
    histo_MC -> GetYaxis() -> SetTitle("");
    
    histo_MC -> SetLineColor(kRed);
    histo_DA -> SetMarkerColor(kBlack);
    histo_DA -> SetMarkerSize(0.5);
    
    histo_MC -> SetFillColor(kRed);
    histo_MC -> SetFillStyle(3003);
        
    histo_MC -> SetLineWidth(2);
    histo_DA -> SetLineWidth(2);
    
    float maximum = -1.;
    if( histo_DA -> GetMaximum() > maximum ) maximum = histo_DA -> GetMaximum();
    if( histo_MC -> GetMaximum() > maximum ) maximum = histo_MC -> GetMaximum();
    histo_DA -> SetMaximum( 1.2*maximum );
    histo_DA -> SetMinimum( 0. );
    histo_MC -> SetMaximum( 1.2*maximum );
    histo_MC -> SetMinimum( 0. );

    TCanvas* c1 = new TCanvas("c1","c1");
    c1 -> cd();
    
    TPaveStats* st_DA = new TPaveStats();
    TPaveStats* st_MC = new TPaveStats();
   
    histo_MC -> Draw("");
    gPad -> Update();
    st_MC= (TPaveStats*)(histo_MC->GetListOfFunctions()->FindObject("stats"));
    st_MC->SetX1NDC(0.82); //new x start position
    st_MC->SetX2NDC(0.99); //new x end position
    st_MC->SetY1NDC(0.56); //new y start position
    st_MC->SetY2NDC(0.68); //new y end position
    st_MC->SetTextColor(kRed);
    st_MC->Draw("sames");

    histo_DA -> Draw("sames");
    gPad -> Update();
    st_DA= (TPaveStats*)(histo_DA->GetListOfFunctions()->FindObject("stats"));
    st_DA->SetX1NDC(0.82); //new x start position
    st_DA->SetX2NDC(0.99); //new x end position
    st_DA->SetY1NDC(0.70); //new y start position
    st_DA->SetY2NDC(0.82); //new y end position
    st_DA->SetTextColor(kBlack);
    st_DA->Draw("sames");
    
    TLegend* legend = new TLegend(0.82, 0.82, 0.99, 0.94);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);

    legend -> AddEntry(histo_DA,label_DA.c_str(),"L");
    legend -> AddEntry(histo_MC,label_MC.c_str(),"F");
    legend -> Draw("same");

    TLatex* latex = new TLatex(0.026, 0.97, xTitle.c_str());
    latex -> SetNDC();
    latex -> SetTextFont(kBlack);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");

    latex = new TLatex(0.76, 0.50, intDA_buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.030);
    latex -> Draw("same");
      
    latex = new TLatex(0.76, 0.44, intMC_buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextColor(kRed);
    latex -> SetTextSize(0.030);
    latex -> Draw("same");

    latex = new TLatex(0.76, 0.38, FoM_buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextColor(kBlack);
    latex -> SetTextSize(0.030);
    latex -> Draw("same");


    c1 -> Print((outputDir+"comp_h_"+xTitle+doCUT+".png").c_str(),"png");
    c1 -> Print((outputDir+"comp_h_"+xTitle+doCUT+".pdf").c_str(),"pdf");

    delete c1;
}

void compareHistosSB(TH1F* histo_DA, TH1F* histo_MC, TH1F* histo_MCbkg,std::string xTitle, std::string label_DA, std::string label_MC,std::string label_MCbkg , std::string outputDir, std::string doCUT, double& integral_data, double& integral_Ref, double& integral, std::string CAT){   
    
    char intDA_buffer[50];
    sprintf(intDA_buffer,("#evt-"+label_DA+" = %1.2e").c_str(),integral_data);

    char intMC_buffer[50];
    sprintf(intMC_buffer,("#evt-"+label_MC+" = %1.2e").c_str(),integral_Ref);

    char intMCbkg_buffer[50];
    sprintf(intMCbkg_buffer,("#evt-"+label_MCbkg+" = %1.2e").c_str(),integral);
    
    histo_DA -> GetXaxis() -> SetLabelSize(0);
    histo_DA -> GetYaxis() -> SetLabelSize(0.04);
    histo_DA -> GetXaxis() -> SetTitleSize(0.05);
    histo_DA -> GetYaxis() -> SetTitleSize(0.05);
    histo_DA -> GetYaxis() -> SetTitleOffset(1.25);

    histo_DA -> GetXaxis() -> SetTitle(xTitle.c_str());
    histo_MC -> GetXaxis() -> SetTitle(xTitle.c_str());
    histo_MCbkg -> GetXaxis() -> SetTitle(xTitle.c_str());

    histo_DA -> GetYaxis() -> SetTitle("");
    histo_MC -> GetYaxis() -> SetTitle("");
    histo_MCbkg -> GetYaxis() -> SetTitle("");
    
    histo_MCbkg -> SetLineColor(kRed);
    histo_MC -> SetLineColor(kBlue);
    histo_DA -> SetMarkerColor(kBlack);
    histo_DA -> SetMarkerSize(0.5);
    
    histo_MCbkg -> SetFillColor(kRed);
    histo_MCbkg -> SetFillStyle(3003);
        
    histo_MCbkg -> SetLineWidth(2);
    histo_DA -> SetLineWidth(2);
   
    float maximum = -1.;
    if( histo_DA -> GetMaximum() > maximum ) maximum = histo_DA -> GetMaximum();
    if( histo_MC -> GetMaximum() > maximum ) maximum = histo_MC -> GetMaximum();
    if( histo_MCbkg -> GetMaximum() > maximum ) maximum = histo_MCbkg -> GetMaximum();
    histo_DA -> SetMaximum( 1.2*maximum );
    histo_DA -> SetMinimum( 0. );
    histo_MC -> SetMaximum( 1.2*maximum );
    histo_MC -> SetMinimum( 0. );
    histo_MCbkg -> SetMaximum( 1.2*maximum );
    histo_MCbkg -> SetMinimum( 0. );
      
    TCanvas* c1 = new TCanvas("c1","c1");
    c1 -> cd();
    
    histo_MCbkg -> Draw(""); 
    histo_MC -> Draw("sames");
    histo_DA -> Draw("e,sames");
    
    TLegend* legend = new TLegend(0.82, 0.82, 0.99, 0.94);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(histo_DA,label_DA.c_str(),"PL");
    legend -> AddEntry(histo_MCbkg,label_MCbkg.c_str(),"F");
    legend -> AddEntry(histo_MC,label_MC.c_str(),"L");
    legend -> Draw("same");
    
    TLatex* latex = new TLatex(0.026, 0.97, xTitle.c_str());
    latex -> SetNDC();
    latex -> SetTextFont(kBlack);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");

    latex = new TLatex(0.76, 0.50, intDA_buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.030);
    latex -> Draw("same");
      
    latex = new TLatex(0.76, 0.44, intMCbkg_buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextColor(kRed);
    latex -> SetTextSize(0.030);
    latex -> Draw("same");

    latex = new TLatex(0.76, 0.38, intMC_buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextColor(kBlue);
    latex -> SetTextSize(0.030);
    latex -> Draw("same");
    
    c1 -> Print((outputDir+"SB_comp_h_"+xTitle+doCUT+CAT+".png").c_str(),"png");
    c1 -> Print((outputDir+"SB_comp_h_"+xTitle+doCUT+CAT+".pdf").c_str(),"pdf");

    delete c1;
}


void draw2DHisto(TH2F* h1, std::string xTitle, std::string yTitle, std::string label, std::string outputDir, std::string doCUT){
	    
    TCanvas* c1 = new TCanvas("c1","c1");
    c1 -> cd();
    
    h1 -> GetXaxis() -> SetLabelSize(0.04);
    h1 -> GetYaxis() -> SetLabelSize(0.04);
    h1 -> GetXaxis() -> SetTitleSize(0.05);
    h1 -> GetYaxis() -> SetTitleSize(0.05);
    h1 -> GetYaxis() -> SetTitleOffset(1.25);

    h1 -> GetXaxis() -> SetTitle(xTitle.c_str());
    h1 -> GetYaxis() -> SetTitle(yTitle.c_str());
    
    h1 -> Draw("colz");
    
    c1 -> Print((outputDir+"comp_h2_"+xTitle+"_vs_"+yTitle+doCUT+"_"+label+".png").c_str(),"png");
    c1 -> Print((outputDir+"comp_h2_"+xTitle+"_vs_"+yTitle+doCUT+"_"+label+".pdf").c_str(),"pdf");
    
    delete c1;
}


void drawNumEvents(TH1D* h1,std::string yTitle,bool& isLog, std::string label,std::string Name,std::string outputDir,std::map<int,std::string> cutString){
	    
    std::string doCUT = cutString[cutString.size()-1];
    
    TCanvas* c1 = new TCanvas("c1","c1");
    c1 -> cd();
    c1 -> SetGridx();
    c1 -> SetGridy();
    if(isLog)c1 -> SetLogy();
    
    h1 -> GetXaxis() -> SetLabelSize(0.02);
    h1 -> GetYaxis() -> SetLabelSize(0.04);
    h1 -> GetXaxis() -> SetTitleSize(0.05);
    h1 -> GetYaxis() -> SetTitleSize(0.05);
    h1 -> GetYaxis() -> SetTitleOffset(1.25);

    //h1 -> GetXaxis() -> SetTitle("cuts");
    h1 -> GetYaxis() -> SetTitle(yTitle.c_str());

    for(unsigned int ii = 0; ii < cutString.size(); ii++){
        std::string cutName = cutString[ii].erase(0,6);
        if(ii == 0)  h1->GetXaxis()->SetBinLabel(ii+1,"noCut");
        else h1->GetXaxis()->SetBinLabel(ii+1,cutName.c_str());
    }
    
    h1 -> SetMarkerStyle(22);
    h1 -> SetMarkerColor(kBlack);
    h1 -> SetMarkerSize(2.);
    h1 -> Draw("P");
    
    c1 -> Print((outputDir+Name+doCUT+"_"+label+".png").c_str(),"png");
    c1 -> Print((outputDir+Name+doCUT+"_"+label+".pdf").c_str(),"pdf");
    
    delete c1;
}

