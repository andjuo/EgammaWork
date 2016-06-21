#include <TROOT.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TLegend.h>
#include <TRandom3.h>
#include "../RooUnfoldInterface/src/RooUnfoldResponse.h"
#include "../RooUnfoldInterface/src/RooUnfoldBayes.h"
#include "../RooUnfoldInterface/src/RooUnfoldInvert.h"
#include <iostream>
//#include "unfold_RMatrix.h"

#include <TMatrixD.h>
#include <TVectorD.h>

// ---------------------------------------------------------------------

void prepareHist(TH1D* h, int color, int markerStyle)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerStyle);
}


// ---------------------------------------------------------------------

void printHistRatio(const TH1D* h1a, const TH1D* h1b, int extraRange=0)
{
  std::cout << "histRatio :\n";
  std::cout << " h1a: " << h1a->GetName() << " " << h1a->GetTitle() << "\n";
  std::cout << " h1b: " << h1b->GetName() << " " << h1b->GetTitle() << "\n";
  int d=(extraRange) ? 1:0;

  for (int ibin=1-d; ibin<=h1a->GetNbinsX()+d; ibin++) {
    if ( (h1a->GetBinLowEdge(ibin) != h1b->GetBinLowEdge(ibin)) ||
	 (h1a->GetBinWidth(ibin) != h1b->GetBinWidth(ibin)) ) {
      std::cout << "bining mismatch at ibin=" << ibin << "\n";
      return;
    }
    std::cout << "ibin=" << ibin << " " << h1a->GetBinLowEdge(ibin)
	      << " " << (h1a->GetBinLowEdge(ibin)+h1a->GetBinWidth(ibin))
	      << "  " << h1a->GetBinContent(ibin) << " +- "
	      << h1a->GetBinError(ibin)
	      << "  " << h1b->GetBinContent(ibin) << " +- "
	      << h1b->GetBinError(ibin)
	      << "  " << h1a->GetBinContent(ibin)/h1b->GetBinContent(ibin)
	      << "\n";
  }
}

// ---------------------------------------------------------------------

//TH1D *do_unfoldByMInv(const TH1D *hMeas, const TH2D *h2Resp);
void drawResponseMatrix(const TH2D *h2, TString tag);

// ---------------------------------------------------------------------

void unfClosureTest(int fsrUnf=0, int rndToy=0, int flag_drawResponseMatrix=0)
{
  TString fname="dyll_M-50.root";
  fname="dyee_resp.root";
  TString measName= (fsrUnf) ? "genPostFSR4pi_Mass" : "recoPostFSRinAcc_Mass";
  TString trueName= (fsrUnf) ? "genPreFSR4pi_Mass" : "genPostFSRinAcc_Mass";
  TString respName= (fsrUnf) ? "RM_FSR4pi" : "RM_postFSR";
  TString rooRespName= (fsrUnf) ? "rooUnf_fsrResp" : "rooUnf_detResResp";
  TString tag= (fsrUnf) ? "_FSR" : "_Resol";

  // use original distributions for the resolution unfolding
  if (fsrUnf==-1) {
    measName="ZMass";
    trueName="genPostFSR_Mass";
    respName="RM_postFSR";
    tag="_origResol";
  }

  // reposnse matrix constructed from scratch
  RooUnfoldResponse *rooResp=NULL;

  // Get histograms
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return;
  }
  TH1D *hMeas= (TH1D*) fin.Get(measName);
  TH1D *hTrue= (TH1D*) fin.Get(trueName);
  TH2D *h2Resp= NULL; //(TH2D*) fin.Get(respName);
  rooResp= (RooUnfoldResponse*) fin.Get(rooRespName);
  hMeas->SetDirectory(0); hTrue->SetDirectory(0);
  if (h2Resp) h2Resp->SetDirectory(0);
  fin.Close();

  rooResp->UseOverflow(false);
  if (!h2Resp) h2Resp= (TH2D*)rooResp->Hresponse();

  std::cout << "number of entries: " << hMeas->Integral() << ", " << hTrue->Integral() << ", " << h2Resp->Integral() << "\n";

  TH1D* hMeas_Mdf=(TH1D*) hMeas->Clone();
  hMeas_Mdf->SetDirectory(0);
  if (rndToy) {
    for (int ibin=1; ibin<=hMeas_Mdf->GetNbinsX(); ibin++) {
      double v= gRandom->Gaus(hMeas_Mdf->GetBinContent(ibin),
			      hMeas_Mdf->GetBinError(ibin));
      hMeas_Mdf->SetBinContent(ibin, v);
      hMeas_Mdf->SetBinError(ibin, 0.);
    }
  }

  // Draw the migration (response) matrix, if needed
  if (flag_drawResponseMatrix) {
    drawResponseMatrix(h2Resp,tag);
    return;
  }

  // Create RooUnfoldBayes for closure test
  int nIters=1;
  TString rubName="RooUnfBayes_" + tag + Form("_%d",nIters);
  RooUnfoldBayes ruBayes( rooResp,hMeas_Mdf,nIters,false,rubName,rubName);
  TH1D* hUnf=(TH1D*)ruBayes.Hreco();
  hUnf->SetDirectory(0);

  RooUnfoldInvert rooMInvert( rooResp,hMeas_Mdf, "rooMInvert","rooMInvert");
  TH1D* hUnfByMInv=(TH1D*)rooMInvert.Hreco();
  hUnfByMInv->SetDirectory(0);


  // Define labels for the legend
  TString labelMeas="measured";
  TString labelTrue="true";
  TString labelUnf =TString("unf.by Bayes") + Form("[%d]",nIters);
  TString labelUnfMInv="unf.by M inv.";

  // rename x label
  //hMeas->GetXaxis()->SetTitle("M_{ee} [GeV]");

  // Make the plot
  TString name="cx_" + tag;
  TCanvas *cx=new TCanvas(name,name,600,600);
  cx->SetLogx();
  cx->SetLogy();
  //hMeas->GetYaxis()->SetRangeUser(1e4,1e7);
  prepareHist(hMeas,kBlack,24);
  prepareHist(hTrue,38,20);
  prepareHist(hUnf,kRed,4); hUnf->SetLineStyle(2);
  hMeas->Draw("LPE");
  hTrue->Draw("LPE same");
  hUnf->Draw("hist same");

  if (hUnfByMInv) {
    prepareHist(hUnfByMInv,kGreen+1,2);
    hUnfByMInv->Draw("LPE same");
  }
  cx->Update();

  TLegend *leg= new TLegend(0.15,0.23,0.35,0.38);
  leg->AddEntry(hMeas,labelMeas,"LP");
  leg->AddEntry(hTrue,labelTrue,"LP");
  leg->AddEntry(hUnf, labelUnf, "LP");
  if (hUnfByMInv) {
    leg->AddEntry(hUnfByMInv, labelUnfMInv, "LP");
  }
  leg->Draw();
  cx->Modified();
  cx->Update();

  if (1) {
    printHistRatio(hTrue, hUnf);
  }

}

// ----------------------------------------------------------
// ----------------------------------------------------------

void drawResponseMatrix(const TH2D* h2_inp, TString tag)
{
  TString h2Name="h2_drawResponseMatrix" + tag;
  TH2D *h2=(TH2D*)h2_inp->Clone(h2Name);
  h2->SetStats(0);

  TString cname="cResponse" + tag;
  TCanvas *cr= new TCanvas(cname,cname,600,600);
  cr->SetRightMargin(0.18);
  cr->SetLogx(1);
  cr->SetLogy(1);
  h2->Draw("COLZ");
  cr->Update();
  //h2->GetZaxis()->SetRangeUser(1,maxVal);
  //cr->Modified();
  //cr->Update();
  return;
}

// ----------------------------------------------------------
// ----------------------------------------------------------
