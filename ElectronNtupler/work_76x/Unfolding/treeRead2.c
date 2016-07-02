#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include "../RooUnfoldInterface/src/RooUnfoldResponse.h"

using std::cout;

//********************************* Unfolding Detector Resolution*********************

const int nBins=43;  
Double_t xbin[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

int insideMassRange(double m)
{
  return ((m<xbin[0]) || (m>xbin[nBins])) ? 0 : 1;
}

void treeRead2()
{
  const int nFiles= 11;
  const TString fnames[nFiles] = {
    "outMass_10to50.root",
    "outMass_50to100.root",
    "outMass_100to200.root",
    "outMass_200to400.root",
    "outMass_400to500.root",
    "outMass_500to700.root",
    "outMass_700to800.root",
    "outMass_800to1000.root",
    "outMass_1500to2000.root",
    "outMass_1000to1500.root",
    "outMass_2000to3000.root"
  };
  TString outFileName="detRes.root";

  double massReco, massGen, lumiWeight, genWeight, PUWeight;
  double massPreFSR;
  char isReco, isGen, BB, BE, EE;

  const int nBins=43;
  Double_t xbin[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
  //Double_t xbin[32] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};
  int rangeInt=int(xbin[nBins]+0.1);

  int i=0;
  TH1D *reco_EEMass[2]; TH1D *gen_EEMass[2]; TH2D *RespMatrix[2]; TH2D *isMiss_RespMatrix[2]; TH2D *isFake_RespMatrix[2];
  TH2D * resp_DetRes[2];
  TH1D * h1_preFSRinAcc_1GeV_wPU;
  TH1D * h1_preFSRinAcc_1GeV_noPU;

  reco_EEMass[i] = new TH1D("reco_EEMass", "reco_EEMass", 43.0, xbin);
  gen_EEMass[i]  = new TH1D("gen_EEMass", "gen_EEMass", 43.0, xbin);
  RespMatrix[i]  = new TH2D("RespMatrix", "RespMatrix", 43.0, xbin, 43.0, xbin);
  isMiss_RespMatrix[i] = new TH2D("isMiss_RespMatrix", "isMass_RespMatrix", 43.0, xbin, 43.0, xbin);
  isFake_RespMatrix[i] = new TH2D("isFake_RespMatrix", "isFake_RespMatrix", 43.0, xbin, 43.0, xbin);

  resp_DetRes[i] = new TH2D("resp_DetRes", "", 43.0, xbin, 43.0, xbin);
  reco_EEMass[i]->Sumw2(); gen_EEMass[i]->Sumw2(); RespMatrix[i]->Sumw2(); isMiss_RespMatrix[i]->Sumw2(); isFake_RespMatrix[i]->Sumw2();

  h1_preFSRinAcc_1GeV_noPU= new TH1D("h1_preFSRinAcc_1GeV_noPU", "preFSR in Acc;M_{ee,preFSR} [GeV];weighed count", rangeInt, 1., rangeInt);
  h1_preFSRinAcc_1GeV_wPU= new TH1D("h1_preFSRinAcc_1GeV_wPU", "preFSR in Acc;M_{ee,preFSR} [GeV];weighed count (wPU)", rangeInt, 1., rangeInt);
  h1_preFSRinAcc_1GeV_noPU->Sumw2(); h1_preFSRinAcc_1GeV_wPU->Sumw2();

  // prototypes to define binning for RooUnfoldResponse objects
  TH1D *h1_reco_EEMass_= new TH1D("h1_reco_EEMass_", "reco_EEMass_;M_{ee,RECO} [GeV];weighed count (wPU)", nBins, xbin);
  TH1D *h1_postFSRinAcc_EEMass_= new TH1D("h1_postFSRinAcc_EEMass_", "postFSR_inAcc_EEMass_;M_{ee,postFSR,inAcc} [GeV];weighted count (wPU)", nBins, xbin);
  TH2D *h2_migMatrixDetRes_= new TH2D("h2_migMatrixDetRes_", "migMatrix for DetRes;M_{ee,RECO} [GeV];M_{ee,postFSR,inAcc} [GeV]",nBins,xbin,nBins,xbin);
  h1_reco_EEMass_->Sumw2(); h1_postFSRinAcc_EEMass_->Sumw2(); h2_migMatrixDetRes_->Sumw2();
  RooUnfoldResponse detRes(h1_reco_EEMass_, h1_postFSRinAcc_EEMass_, h2_migMatrixDetRes_, "detRes" );

  // stand-alone distributions
  TH1D* h1_reco_EEMass= (TH1D*)h1_reco_EEMass_->Clone("h1_reco_EEMass");
  TH1D* h1_postFSRinAcc_EEMass= (TH1D*)h1_postFSRinAcc_EEMass_->Clone("h1_postFSRinAcc_EEMass");


  for(int iFile=0; iFile<nFiles; iFile++){

    TFile *f = TFile::Open(fnames[iFile]);
    if (!f->IsOpen()) {
      std::cout << "failed to open the file <" << fnames[iFile] << ">\n";
      return;
    }

    TTree *t = (TTree*)f->Get("tree");
    Long64_t entries = t->GetEntries();
    cout<< f->GetName() << " num.entries: "<<entries<<endl;

    t->SetBranchAddress("massReco",&massReco);
    t->SetBranchAddress("massGen",&massGen);
    t->SetBranchAddress("massPreFSR",&massPreFSR);
    t->SetBranchAddress("lumiWeight",&lumiWeight);
    t->SetBranchAddress("genWeight",&genWeight);
    t->SetBranchAddress("PUWeight",&PUWeight);
    t->SetBranchAddress("isReco",&isReco);
    t->SetBranchAddress("isGen",&isGen);
    t->SetBranchAddress("BB",&BB);
    t->SetBranchAddress("BE",&BE);
    t->SetBranchAddress("EE",&EE);

    for(Long64_t j=0;j<entries;j++){
      t->GetEntry(j);

      double scale = lumiWeight*genWeight*PUWeight*2316.969;
      double scale_noPU= scale/PUWeight;

      // control the distributions
      h1_preFSRinAcc_1GeV_noPU->Fill(massPreFSR,scale_noPU);
      h1_preFSRinAcc_1GeV_wPU->Fill(massPreFSR,scale);


      // prepare the response matrix by RCh
      if(massReco > -999.) reco_EEMass[i]->Fill(massReco,scale);
      if(massGen > -999.) gen_EEMass[i]->Fill(massGen,scale);

      if(isReco && isGen) RespMatrix[i]->Fill(massReco,massGen,scale);
      if(!isReco && isGen) isMiss_RespMatrix[i]->Fill(massReco,massGen,scale);
      if(isReco && !isGen) isFake_RespMatrix[i]->Fill(massReco,massGen,scale);

      if(massReco > -999. && massGen > -999.) {
      	if (insideMassRange(massReco) || insideMassRange(massGen)) {
      		//if (!insideMassRange(massReco) &&  insideMassRange(massGen)) resp_DetRes->Fill(massReco, massGen, scale);
      		//else if ( insideMassRange(MassReco) && !insideMassRange(massGen)) resp_DetRes->Fill(massReco, massGen, scale);

      		resp_DetRes[i]->Fill(massReco, massGen, scale);
        }
      }

      // prepare the response matrix by AJ
      if (massReco > -999.) h1_reco_EEMass->Fill(massReco, scale);
      if (massGen > -999.) h1_postFSRinAcc_EEMass->Fill(massGen, scale);

      if (insideMassRange(massReco) || insideMassRange(massGen)) {
	if (!insideMassRange(massReco)) detRes.Miss(massGen, scale);
	else if (!insideMassRange(massGen)) detRes.Fake(massReco, scale);
	else detRes.Fill(massReco,massGen, scale);
      }

    }

    delete f;
  }


  // Save histograms
  TString outFName="dyee_unf_input.root";
  TFile file2(outFName,"recreate");
  h1_preFSRinAcc_1GeV_wPU->Write();
  h1_preFSRinAcc_1GeV_noPU->Write();
  h1_reco_EEMass->Write();
  h1_postFSRinAcc_EEMass->Write();
  detRes.Write();
  reco_EEMass[i]->Write();
  gen_EEMass[i]->Write();
  RespMatrix[i]->Write();
  resp_DetRes[i]->Write();
  isMiss_RespMatrix[i]->Write();
  isFake_RespMatrix[i]->Write();
  file2.Write();
  file2.Close();

}
