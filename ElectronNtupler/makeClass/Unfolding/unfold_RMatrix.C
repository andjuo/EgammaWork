#define unfold_RMatrix_cxx
#include "unfold_RMatrix.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>

<<<<<<< Updated upstream
#define TEST_RooUnfoldResponse
=======
//#define TEST_RooUnfoldResponse
>>>>>>> Stashed changes

#ifdef TEST_RooUnfoldResponse
#include "../RooUnfoldInterface/src/RooUnfoldResponse.h"
#endif


Double_t unfold_RMatrix::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

Double_t unfold_RMatrix::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t unfold_RMatrix::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

//bool unfold_RMatrix::mySortFnx(Double_t i, Double_t j) { return i>j; }

void unfold_RMatrix::Loop()
{

#ifdef TEST_RooUnfoldResponse
  RooUnfoldResponse ruResp(nBinsFlat,-0.5,nBinsFlat-0.5,"RUResponse","RUResponse");
#endif


  TFile *file = new TFile("dyll_M-50.root", "recreate");

  int good_elec, gen_elec_pre, gen_elec_post;
  double dR;
  double sum_weights;
  double Z_Mass, Z_Y, Z_Rap, Z_Eta, Z_Pt, Z_Phi, ZMass_G_post, ZMass_R_post, ZMass_di_gPost;
  double ZMass_R_post_barrel, ZMass_R_post_endcap, ZMass_R_post_b_e;

  TLorentzVector ele1,ele2,dielectron;
  TLorentzVector reco1_post,reco2_post,direco_post;

  TLorentzVector reco1_post_barrel,reco2_post_barrel,direco_post_barrel;
  TLorentzVector reco1_post_endcap,reco2_post_endcap,direco_post_endcap;
  TLorentzVector reco1_post_b_e,reco2_post_b_e,direco_post_b_e;

  TLorentzVector gen1_post,gen2_post,digen_post;
  TLorentzVector gPost1,gPost2,di_gPost;

  vector <double> newelePt; vector <double> neweleEta; vector <double> neweleEnr; vector <double> newelePhi; vector <double> neweleCharge;

  vector <double> newgPost_Pt; vector <double> newgPost_Eta; vector <double> newgPost_Enr; vector <double> newgPost_Phi;

  vector <double> recoPost_Pt; vector <double> recoPost_Eta; vector <double> recoPost_Enr; vector <double> recoPost_Phi;
  vector <double> genPost_Pt; vector <double> genPost_Eta; vector <double> genPost_Enr; vector <double> genPost_Phi;

  Double_t xbins[48] = {0,10,15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1200,1500,2000,3000};

  TH1D *genPostFSR_Mass = new TH1D("genPostFSR_Mass", "genPostFSR_Mass;M_{ee,gen,postFSR,inAcc} [GeV];weighted count", 47, xbins);
  TH2D *responsePost = new TH2D("RM_postFSR", "RM_postFSR;M_{ee,gen,inAcc} [GeV];M_{ee,reco,inAcc} [GeV];weighted count", 47, xbins, 47, xbins);
  TH1D *reco_ZMass = new TH1D("ZMass", "ZMass;M_{ee,reco,inAcc} [GeV];weighted count", 47, xbins);

  TH1D *recoPostFSRinAcc_Mass = new TH1D("recoPostFSRinAcc_Mass","recoPostFSRinAcc_Mass;M_{postFSR,reco} [GeV];weighted count", 47, xbins);
  TH1D *genPostFSRinAcc_Mass = new TH1D("genPostFSRinAcc_Mass","genPostFSRinAcc_Mass;M_{postFSR,gen} [GeV];weighted count", 47,xbins);

  TH1D *genPostFSR4pi_Mass= new TH1D("genPostFSR4pi_Mass", "getPostFSR4pi_Mass;M_{postFSR,fullSp} [GeV];weighted count", 47, xbins);
  TH1D *genPreFSR4pi_Mass = new TH1D("genPreFSR4pi_Mass", "genPreFSR4pi_Mass;Mass_{preFSR,fullSp} [GeV];weighted count", 47, xbins);
  TH2D *fsrResponse4pi= new TH2D("RM_FSR4pi","RM_FSR4pi;M_{ee,gen,preFSR4pi} [GeV];M_{ee,gen,postFSR4pi} [GeV];weighted count", 47, xbins, 47, xbins);

  responsePost->Sumw2(); genPostFSR_Mass->Sumw2(); reco_ZMass->Sumw2();
  genPostFSR4pi_Mass->Sumw2(); genPreFSR4pi_Mass->Sumw2(); fsrResponse4pi->Sumw2();

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 500;
  cout<<"entries: "<<nentries<<endl;

  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    // Sorting Reco level
    int index[pt->size()];
    float ptnew[pt->size()];

    for(unsigned int ele=0; ele<pt->size(); ele++)
    {
      ptnew[ele]=pt->at(ele);
    }

    int size = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(size,ptnew,index,true);

    // Sorting Post Gen level
    int index2[gen_postFSR_pt->size()];
    float pt2[gen_postFSR_pt->size()];

    for(unsigned int b=0; b<gen_postFSR_pt->size(); b++)
    {
      pt2[b]=gen_postFSR_pt->at(b);
    }

    int sizen = sizeof(pt2)/sizeof(pt2[0]);
    TMath::Sort(sizen,pt2,index2,true);

    Z_Mass=0.; Z_Y=0.; Z_Rap=0.; Z_Eta=0.; Z_Pt=0.; Z_Phi=0.;
    ZMass_G_post=0.; ZMass_R_post=0.; ZMass_di_gPost=0.;

    good_elec = 0;
    gen_elec_post = 0;
    dR = 0.;

    unsigned int matchGen1 = 0;
    unsigned int matchGen2 = 0;
    unsigned int matchReco1 = 0;
    unsigned int matchReco2 = 0;

    newelePt.clear(); neweleEta.clear(); neweleEnr.clear(); newelePhi.clear(); neweleCharge.clear();

    newgPost_Pt.clear(); newgPost_Eta.clear(); newgPost_Enr.clear(); newgPost_Phi.clear();

    recoPost_Pt.clear(); recoPost_Eta.clear(); recoPost_Enr.clear(); recoPost_Phi.clear(); 
    genPost_Pt.clear(); genPost_Eta.clear(); genPost_Enr.clear(); genPost_Phi.clear();

    //cout<<"1"<<endl;

    if(jentry==26194585) continue;
    
    if(!tauFlag){
      if(nEle>=2) {
	if(single_Ele27) {

	  for(int k=0;k<nEle;k++){

	    if(passMediumId->at(index[k]) == 1){ // && eleEcalDrivenSeed->at(index[k]) == 1)
	      if(fabs(eta->at(index[k])) < 2.5 && !(fabs(etaSC->at(index[k])) > 1.4442 && fabs(etaSC->at(index[k])) < 1.566)){

		if(passMediumId->at(index[k]) == 0) cout<<"Wrong ID: "<<endl;
		//if(eleEcalDrivenSeed->at(index[k]) == 0) cout<<"Not ECAL Driven: "<<endl;

		good_elec = good_elec + 1;

		newelePt.push_back(pt->at(index[k]));
		neweleEta.push_back(eta->at(index[k]));
		neweleEnr.push_back(energy->at(index[k]));
		newelePhi.push_back(phi->at(index[k]));
		neweleCharge.push_back(charge->at(index[k]));

	      } //kin
	    } // ID
	  } // k loop
	} //trigger
      } // nEle
    } // tau flag

    //cout<<"2"<<endl;

    if(good_elec==2){
      if(newelePt.at(0) < newelePt.at(1)) cout<<"event: "<<jentry<<"   "<<"Sorting not proper: "<<"   "<<"reco pt lead: "<<newelePt.at(0)<<"   "<<"reco pt sublead: "<<newelePt.at(1)<<endl;

      if(newelePt.at(0) > 30. && newelePt.at(1) > 10.){

	ele1.SetPtEtaPhiE(newelePt.at(0),neweleEta.at(0),newelePhi.at(0),neweleEnr.at(0));
	ele2.SetPtEtaPhiE(newelePt.at(1),neweleEta.at(1),newelePhi.at(1),neweleEnr.at(1));

	dielectron=ele1+ele2;
	Z_Mass = dielectron.M();

	reco_ZMass->Fill(Z_Mass,theWeight);
	//reco_ZMass->SetBinContent(45,1);

      }
    }

    //cout<<"3"<<endl;

    if(gen_postFSR_pt->size()>=2.){   
      for(unsigned int m=0;m<gen_postFSR_eta->size();m++){

	if(fabs(gen_postFSR_eta->at(index2[m])) < 2.5 && !(fabs(gen_postFSR_eta->at(index2[m])) > 1.4442 && fabs(gen_postFSR_eta->at(index2[m])) < 1.566)){

	  gen_elec_post = gen_elec_post + 1;

	  newgPost_Pt.push_back(gen_postFSR_pt->at(index2[m]));
	  newgPost_Eta.push_back(gen_postFSR_eta->at(index2[m]));
	  newgPost_Enr.push_back(gen_postFSR_ene->at(index2[m]));
	  newgPost_Phi.push_back(gen_postFSR_phi->at(index2[m]));
	}
      }
    }

    if(gen_elec_post==2){

      if(newgPost_Pt.at(0) < newgPost_Pt.at(1)) cout<<"entry: "<<jentry<<"   "<<"Sorting not done properly"<<"   "<<"new gen post pt lead: "<<newgPost_Pt.at(0)<<"   "<<"new gen post pt sublead: "<<newgPost_Pt.at(1)<<endl;

      if(newgPost_Pt.at(0) > 30. && newgPost_Pt.at(1) > 10.){

	gPost1.SetPtEtaPhiE(newgPost_Pt.at(0),newgPost_Eta.at(0),newgPost_Phi.at(0),newgPost_Enr.at(0));
	gPost2.SetPtEtaPhiE(newgPost_Pt.at(1),newgPost_Eta.at(1),newgPost_Phi.at(1),newgPost_Enr.at(1));

	di_gPost=gPost1+gPost2;
	ZMass_di_gPost=di_gPost.M();
	genPostFSR_Mass->Fill(ZMass_di_gPost,theWeight);
	//genPostFSR_Mass->SetBinContent(45,1);

      }
    }

    // Reco-gen matching
    if(gen_elec_post==2){
      for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
	double dR_comp_post = 1000.;
	for(unsigned int igen = 0; igen < newgPost_Eta.size(); igen++){
	  dR = deltaR(neweleEta[ireco], newelePhi[ireco], newgPost_Eta[igen], newgPost_Phi[igen]);

	  if(dR < 0.1  ){
	    if (dR < dR_comp_post)
	    {
	      dR_comp_post = dR;
	      matchGen2 = igen ; matchReco2 = ireco ; //iter[ireco]=igen;
	    }

	    genPost_Pt.push_back(newgPost_Pt[matchGen2]); genPost_Eta.push_back(newgPost_Eta[matchGen2]); genPost_Enr.push_back(newgPost_Enr[matchGen2]); genPost_Phi.push_back(newgPost_Phi[matchGen2]);
	    recoPost_Pt.push_back(newelePt[matchReco2]); recoPost_Eta.push_back(neweleEta[matchReco2]); recoPost_Enr.push_back(neweleEnr[matchReco2]); recoPost_Phi.push_back(newelePhi[matchReco2]);
	  }
	}
      }
    }

    // recalculate values
    ZMass_R_post=0; ZMass_G_post=0;

    if(recoPost_Pt.size()==2)
    {
      reco1_post.SetPtEtaPhiE(recoPost_Pt.at(0),recoPost_Eta.at(0),recoPost_Phi.at(0),recoPost_Enr.at(0));
      reco2_post.SetPtEtaPhiE(recoPost_Pt.at(1),recoPost_Eta.at(1),recoPost_Phi.at(1),recoPost_Enr.at(1));

      direco_post=reco1_post+reco2_post;
      ZMass_R_post = direco_post.M();

    }

    if(genPost_Pt.size()==2)
    {
      gen1_post.SetPtEtaPhiE(genPost_Pt.at(0),genPost_Eta.at(0),genPost_Phi.at(0),genPost_Enr.at(0));
      gen2_post.SetPtEtaPhiE(genPost_Pt.at(1),genPost_Eta.at(1),genPost_Phi.at(1),genPost_Enr.at(1));
      digen_post=gen1_post+gen2_post;
      ZMass_G_post = digen_post.M();
    }

    recoPostFSRinAcc_Mass->Fill(ZMass_R_post, theWeight);
    genPostFSRinAcc_Mass->Fill(ZMass_G_post, theWeight);
    responsePost->Fill(ZMass_G_post,ZMass_R_post,theWeight);

#ifdef TEST_RooUnfoldResponse
    int recoFI= flat_index(ZMass_R_post);
    int genFI = flat_index(ZMass_G_post);
    if ((genFI!=-1) && (recoFI!=-1)) ruResp.Fill(recoFI,genFI,theWeight);
    else if ((genFI==-1) && (recoFI!=-1)) ruResp.Fake(recoFI, theWeight);
    else if ((genFI!=-1) && (recoFI==-1)) ruResp.Miss(genFI, theWeight);
#endif


    // FSR response
    if (!tauFlag) {
      if (nEle>=2) {

	ele1.SetPtEtaPhiE(pt->at(0),eta->at(0),phi->at(0),energy->at(0));
	ele2.SetPtEtaPhiE(pt->at(1),eta->at(1),phi->at(1),energy->at(1));
	dielectron=ele1+ele2;
	Double_t Z_Mass_postFSR4pi = dielectron.M();
	Double_t Z_Mass_preFSR4pi = 0;

	int nPreFsr=gen_preFSR->GetEntriesFast();
	//std::cout << "here are nPreFsr=" << nPreFsr << " electrons\n";
	if (nPreFsr>=2) {
	  int gen_index[nPreFsr];
	  float gen_E[nPreFsr];
	  for (int i=0; i<nPreFsr; i++) {
	    const TLorentzVector *e= (TLorentzVector*)(*gen_preFSR)[i];
	    std::cout << " i = " << e->E() << "\n";
	    gen_E[i]= e->E();
	  }
	  TMath::Sort(int(sizeof(gen_E)/sizeof(gen_E[0])), gen_E,gen_index,true);

	  // Assume that 2 leading gen electrons are of our interest
	  const TLorentzVector *gen_e1= (TLorentzVector*)(*gen_preFSR)[gen_index[0]];
	  const TLorentzVector *gen_e2= (TLorentzVector*)(*gen_preFSR)[gen_index[1]];
	  std::cout << "leading electrons pT: " << gen_e1->Perp() << ", " << gen_e2->Perp() << "\n";
	  std::cout << "leading electrons E: " << gen_e1->E() << ", " << gen_e2->E() << "\n";
	  dielectron=(*gen_e1) + (*gen_e2);
	  Z_Mass_preFSR4pi= dielectron.M();
	}

	genPostFSR4pi_Mass->Fill(Z_Mass_postFSR4pi, theWeight);
	genPreFSR4pi_Mass->Fill (Z_Mass_preFSR4pi, theWeight);
	fsrResponse4pi->Fill(Z_Mass_preFSR4pi,Z_Mass_postFSR4pi,theWeight);
      }
    }


  } // event

#ifdef TEST_RooUnfoldResponse
  ruResp.Write();
#endif

  file->Write();
  file->Close();
}
