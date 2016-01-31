#include <TROOT.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TLegend.h>
#include "../RooUnfoldInterface/src/RooUnfoldResponse.h"
#include "../RooUnfoldInterface/src/RooUnfoldBayes.h"
#include <iostream>

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

TH1D *do_unfoldByMInv(const TH1D *hMeas, const TH2D *h2Resp);
void drawResponseMatrix(const TH2D *h2, TString tag);

// ---------------------------------------------------------------------

void unfClosureTest(int fsrUnf=0, int flag_drawResponseMatrix=0)
{
  TString fname="dyll_M-50.root";
  TString measName= (fsrUnf) ? "genPostFSR4pi_Mass" : "recoPostFSRinAcc_Mass";
  TString trueName= (fsrUnf) ? "genPreFSR4pi_Mass" : "genPostFSRinAcc_Mass";
  TString respName= (fsrUnf) ? "RM_FSR4pi" : "RM_postFSR";
  TString tag= (fsrUnf) ? "_FSR" : "_Resol";

  // use original distributions for the resolution unfolding
  if (fsrUnf==-1) {
    measName="ZMass";
    trueName="genPostFSR_Mass";
    respName="RM_postFSR";
    tag="_origResol";
  }


  // Get histograms
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return;
  }
  TH1D *hMeas= (TH1D*) fin.Get(measName);
  TH1D *hTrue= (TH1D*) fin.Get(trueName);
  TH2D *h2Resp=(TH2D*) fin.Get(respName);
  hMeas->SetDirectory(0); hTrue->SetDirectory(0); h2Resp->SetDirectory(0);
  fin.Close();

  TH2D* h2RespT=(TH2D*)h2Resp->Clone("h2RespTransposed");
  h2RespT->GetXaxis()->SetTitle(h2Resp->GetYaxis()->GetTitle());
  h2RespT->GetYaxis()->SetTitle(h2Resp->GetXaxis()->GetTitle());
  for (int ibin=1; ibin<=h2Resp->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2Resp->GetNbinsY(); jbin++) {
      h2RespT->SetBinContent(jbin,ibin, h2Resp->GetBinContent(ibin,jbin));
      h2RespT->SetBinError  (jbin,ibin, h2Resp->GetBinError  (ibin,jbin));
    }
  }


  std::cout << "number of entries: " << hMeas->Integral() << ", " << hTrue->Integral() << ", " << h2Resp->Integral() << "\n";

  // Create RooUnfold response matrix
  int useTransposed=1;
  if (useTransposed) std::cout << "\n\n\tusing transposed response matrix\n\n";
  TString ruName="RooUnfoldResp";
  if (useTransposed) ruName.Append("T");
  ruName+=tag;
  TH2D *h2ResponseUsed=(useTransposed) ? h2RespT : h2Resp;

  // Draw the migration (response) matrix, if needed
  if (flag_drawResponseMatrix) {
    drawResponseMatrix(h2ResponseUsed,tag);
    return;
  }


  RooUnfoldResponse ruResp(hMeas,hTrue,h2ResponseUsed, ruName,ruName);

  // Create RooUnfoldBayes for closure test
  int nIters=3;
  TString rubName="RooUnfBayes_" + tag + Form("_%d",nIters);
  RooUnfoldBayes ruBayesCT( &ruResp,hMeas,nIters,false,rubName,rubName);
  TH1D* hUnf=(TH1D*)ruBayesCT.Hreco();
  hUnf->SetDirectory(0);


  TH1D* hUnfByMInv=NULL;
  if (0) {
    hUnfByMInv=do_unfoldByMInv(hMeas,h2ResponseUsed);
    if (!hUnfByMInv) return;
  }

  // eliminate influence of 0-15
  hMeas->SetBinContent(1,1e-2);
  hTrue->SetBinContent(1,1e-2);
  hUnf->SetBinContent(1,1e-2);
  if (hUnfByMInv) hUnfByMInv->SetBinContent(1,1e-2);

  // Define labels for the legend
  TString labelMeas="measured";
  TString labelTrue="true";
  TString labelUnf =TString("unf.by Bayes") + Form("[%d]",nIters);
  TString labelUnfMInv="unf.by M inv.";

  // rename x label
  hMeas->GetXaxis()->SetTitle("M_{ee} [GeV]");

  // Make the plot
  TString name="cx_" + tag;
  TCanvas *cx=new TCanvas(name,name,600,600);
  cx->SetLogx();
  cx->SetLogy();
  hMeas->GetYaxis()->SetRangeUser(1e4,1e7);
  prepareHist(hMeas,kBlack,24);
  prepareHist(hTrue,kBlue,20);
  prepareHist(hUnf,kRed,4);
  hMeas->Draw("LPE");
  hTrue->Draw("LPE same");
  hUnf->Draw("hist same");
  if (hUnfByMInv) {
    prepareHist(hUnfByMInv,kGreen+1,5);
    hUnfByMInv->Draw("LPE same");
  }
  cx->Update();

  TLegend *leg= new TLegend(0.15,0.73,0.35,0.88);
  leg->AddEntry(hMeas,labelMeas,"LP");
  leg->AddEntry(hTrue,labelTrue,"LP");
  leg->AddEntry(hUnf, labelUnf, "LP");
  if (hUnfByMInv) {
    leg->AddEntry(hUnfByMInv, labelUnfMInv, "LP");
  }
  leg->Draw();
  cx->Modified();
  cx->Update();

}

// ----------------------------------------------------------
// ----------------------------------------------------------

void drawResponseMatrix(const TH2D* h2_inp, TString tag)
{
  TString h2Name="h2_drawResponseMatrix" + tag;
  TH2D *h2=(TH2D*)h2_inp->Clone(h2Name);
  h2->SetStats(0);

  // clear lower edges
  if (1) {
    std::cout << "\n\tClearing edges for the plot\n";
    h2->SetBinContent(1,1,1e-4);
    for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
      h2->SetBinContent(1,ibin,1e-4);
      h2->SetBinContent(ibin,1,1e-4);
    }
  }

  double maxVal=0;
  int ibin0=0, jbin0=0;
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      double v= h2->GetBinContent(ibin,jbin);
      if (v > maxVal) {
	if (v>1e7) std::cout << "at " << ibin << "," << jbin << " " << v << "\n";
	ibin0=ibin; jbin0=jbin;
	maxVal=v;
      }
    }
  }
  std::cout << "maxVal=" << maxVal << " at " << ibin0 << "," << jbin0 << "\n";

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

TH1D *do_unfoldByMInv(const TH1D *hMeas, const TH2D *h2Resp)
{
  TMatrixD M(h2Resp->GetNbinsX(), h2Resp->GetNbinsY());
  TVectorD vInp(hMeas->GetNbinsX());

  for (int ibin=1; ibin<=hMeas->GetNbinsX(); ++ibin) {
    vInp[ibin-1] = hMeas->GetBinContent(ibin);
  }

  for (int ibin=1; ibin<=h2Resp->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2Resp->GetNbinsY(); jbin++) {
      M(ibin-1,jbin-1) = h2Resp->GetBinContent(jbin,ibin); // note index exchange!!
      // remove underflow
      if ((ibin==1) || (jbin==1)) M(ibin-1,jbin-1)=0.;
    }
  }

  // check that the matrix is invertable
  for (int i=0; i<M.GetNrows(); i++) {
    double sumR=0;
    double sumC=0;
    for (int j=0; j<M.GetNcols(); j++) {
      sumR += M(i,j)*M(i,j);
      sumC += M(j,i)*M(j,i);
    }
    // check if correction is needed
    if ((sumR==0) && (sumC==0)) {
      M(i,i)=1.;
    }
  }

  // Draw migration
  if (0) {
    TString cname="cm";
    TCanvas *cm=new TCanvas(cname,cname,600,600);
    M.Draw("COLZ");
    cm->Update();
    return NULL;
  }

  // Normalize
  if (0) {
    for (int ir=0; ir<M.GetNrows(); ir++) {
      double sumR=0;
      for (int ic=0; ic<M.GetNcols(); ic++) {
	sumR+= M(ir,ic);
      }
      if (sumR==0) M(ir,ir)=1.;
      else
      for (int ic=0; ic<M.GetNcols(); ic++) {
	M(ir,ic) /= sumR;
      }
    }
  }
  else {
    for (int ic=0; ic<M.GetNcols(); ic++) {
      double sumC=0;
      for (int ir=0; ir<M.GetNrows(); ir++) {
	sumC+= M(ir,ic);
      }
      if (sumC==0) M(ic,ic)=1.;
      else
      for (int ir=0; ir<M.GetNrows(); ir++) {
	M(ir,ic) /= sumC;
      }
    }
  }

  // Draw probabilistic normalize
  if (0) {
    TString cname="cmProb";
    TCanvas *cm=new TCanvas(cname,cname,600,600);
    M.Draw("COLZ");
    cm->Update();
    return NULL;
  }



  // unfold by matrix inversion
  TMatrixD Minv(M);
  Double_t det=0;
  M.Invert(&det);
  std::cout << "det=" << det << "\n";
  TVectorD vOut= Minv*vInp;

  // construct the histogram
  TH1D* hOut=(TH1D*)hMeas->Clone("hUnfByInv");
  hOut->SetDirectory(0);

  for (int ibin=1; ibin<=hOut->GetNbinsX(); ibin++) {
    hOut->SetBinContent(ibin, vOut[ibin-1]);
    std::cout << "ibin=" << ibin << ", vOut=" << vOut[ibin-1] << "\n"; 
    // error determination needs some sophistication for error propagation
    hOut->SetBinError  (ibin, 0.);
  }

  return hOut;
}

// ----------------------------------------------------------
