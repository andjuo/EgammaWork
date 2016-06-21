#include "unfold_RMatrix.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <math.h>
#include "../RooUnfoldInterface/src/RooUnfoldResponse.h"

const int nBins=43;
Double_t xbins[nBins+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
			   64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
			   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
			   200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
			   830, 1000, 1500, 3000 };

inline
int insideMassRange(double m)
{
  return ((m<xbins[0]) || (m>xbins[nBins])) ? 0 : 1;
}


// -----------------------------------------------------------------

void createRespM(Long64_t maxEntries_user=10)
{

  /*
  TFile fin("test/DY_check.root");
  if (!fin.IsOpen()) { std::cout << "failed to open the file\n"; return; }
  TTree *t=(TTree*)fin.Get("ntupler/ElectronTree");
  if (!t) { std::cout << "failed to get the tree\n"; return; }
  std::cout << "got tree" << std::endl;
  unfold_RMatrix R(t);
  if (R.fChain == 0) return;
  */

  TString path="/tmp/andriusj/";
  TString inpFiles=path+"DY_10to50.root "+ path+"DY_50toInf.root";
  unfold_RMatrix R(inpFiles,"ntupler/ElectronTree");
  if (R.fChain==0) return;

  TFile *file = new TFile("dyee_resp.root", "recreate");

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

  TH1D *genPostFSR_Mass = new TH1D("genPostFSR_Mass", "genPostFSR_Mass;M_{ee,gen,postFSR,inAcc} [GeV];weighted count", nBins, xbins);
  TH2D *responsePost = new TH2D("RM_postFSR", "RM_postFSR;M_{ee,gen,inAcc} [GeV];M_{ee,reco,inAcc} [GeV];weighted count", nBins, xbins, nBins, xbins);
  TH1D *reco_ZMass = new TH1D("ZMass", "ZMass;M_{ee,reco,inAcc} [GeV];weighted count", nBins, xbins);

  TH1D *recoPostFSRinAcc_Mass = new TH1D("recoPostFSRinAcc_Mass","recoPostFSRinAcc_Mass;M_{postFSR,reco} [GeV];weighted count", nBins, xbins);
  TH1D *genPostFSRinAcc_Mass = new TH1D("genPostFSRinAcc_Mass","genPostFSRinAcc_Mass;M_{postFSR,gen} [GeV];weighted count", nBins,xbins);

  TH1D *genPostFSR4pi_Mass= new TH1D("genPostFSR4pi_Mass", "getPostFSR4pi_Mass;M_{postFSR,fullSp} [GeV];weighted count", nBins, xbins);
  TH1D *genPreFSR4pi_Mass = new TH1D("genPreFSR4pi_Mass", "genPreFSR4pi_Mass;Mass_{preFSR,fullSp} [GeV];weighted count", nBins, xbins);
  TH2D *fsrResponse4pi= new TH2D("RM_FSR4pi","RM_FSR4pi;M_{ee,gen,preFSR4pi} [GeV];M_{ee,gen,postFSR4pi} [GeV];weighted count", 47, xbins, 47, xbins);

  responsePost->Sumw2(); genPostFSR_Mass->Sumw2(); reco_ZMass->Sumw2();
  genPostFSR4pi_Mass->Sumw2(); genPreFSR4pi_Mass->Sumw2(); fsrResponse4pi->Sumw2();


  // prototypes for binning
  TH1D* h1_genReco_proto=(TH1D*)reco_ZMass->Clone("h1_genReco_proto");
  TH1D* h1_genPostFsrInAcc_proto=(TH1D*)genPostFSR_Mass->Clone("h1_genPostFsrInAcc_proto");
  TH2* h2_genMigMatrix_proto=(TH2D*)responsePost->Clone("h2_genMigMatrix_proto");
  h1_genReco_proto->SetDirectory(0);
  h1_genPostFsrInAcc_proto->SetDirectory(0);
  h2_genMigMatrix_proto->SetDirectory(0);
  RooUnfoldResponse* resp_DetRes= new RooUnfoldResponse(h1_genReco_proto,h1_genPostFsrInAcc_proto,h2_genMigMatrix_proto, "rooUnf_detResResp", "detResResp;M_{ee} [GeV];det.res. unfolded yield");

  TH1D* h1_genPostFSR4pi_Mass_proto= (TH1D*)genPostFSR4pi_Mass->Clone("h1_genPostFSR4pi_Mass_proto");
  TH1D* h1_genPreFSR4pi_Mass_proto=(TH1D*)genPreFSR4pi_Mass->Clone("h1_genPreFSR4pi_Mass_proto");
  TH2D* h2_genFSRMigMatrix_proto=(TH2D*)fsrResponse4pi->Clone("h2_genFSRMigMatrix_proto");
  h1_genPostFSR4pi_Mass_proto->SetDirectory(0);
  h1_genPreFSR4pi_Mass_proto->SetDirectory(0);
  h2_genFSRMigMatrix_proto->SetDirectory(0);
  RooUnfoldResponse* resp_FSR= new RooUnfoldResponse(h1_genPostFSR4pi_Mass_proto,h1_genPreFSR4pi_Mass_proto,h2_genFSRMigMatrix_proto,"rooUnf_fsrResp", "FSRResp;M_{ee} [GeV];FSR unfolded yield");

  // Read data

  Long64_t nentries = R.fChain->GetEntriesFast();
  //Long64_t nentries = 500;
  cout<<"entries: "<<nentries<<endl;
  if ((maxEntries_user>0) && (nentries>maxEntries_user)) {
    nentries=maxEntries_user;
    std::cout << "only " << nentries << " will be processed\n";
  }

  sum_weights = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = R.LoadTree(jentry);
    if (ientry < 0) break;
    nb = R.fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%1000000 == 0){
      cout << "Events Processed :  " << jentry << endl;
    }

    // Sorting Reco level
    int index[R.pt->size()];
    float ptnew[R.pt->size()];

    for(unsigned int ele=0; ele<R.pt->size(); ele++)
    {
      ptnew[ele]=R.pt->at(ele);
    }

    int size = sizeof(ptnew)/sizeof(ptnew[0]);
    TMath::Sort(size,ptnew,index,true);

    // Sorting Post Gen level
    int index2[R.gen_postFSR_pt->size()];
    float pt2[R.gen_postFSR_pt->size()];

    for(unsigned int b=0; b<R.gen_postFSR_pt->size(); b++)
    {
      pt2[b]=R.gen_postFSR_pt->at(b);
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

    //if(jentry==26194585) continue;
    
    if(!R.tauFlag){
      if(R.nEle>=2) {
	if(R.single_Ele27) {

	  for(int k=0;k<R.nEle;k++){

	    if(R.passMediumId->at(index[k]) == 1){ // && eleEcalDrivenSeed->at(index[k]) == 1)
	      if(fabs(R.eta->at(index[k])) < 2.5 && !(fabs(R.etaSC->at(index[k])) > 1.4442 && fabs(R.etaSC->at(index[k])) < 1.566)){

		if(R.passMediumId->at(index[k]) == 0) cout<<"Wrong ID: "<<endl;
		//if(eleEcalDrivenSeed->at(index[k]) == 0) cout<<"Not ECAL Driven: "<<endl;

		good_elec = good_elec + 1;

		newelePt.push_back(R.pt->at(index[k]));
		neweleEta.push_back(R.eta->at(index[k]));
		neweleEnr.push_back(R.energy->at(index[k]));
		newelePhi.push_back(R.phi->at(index[k]));
		neweleCharge.push_back(R.charge->at(index[k]));

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

	reco_ZMass->Fill(Z_Mass,R.theWeight);
	//reco_ZMass->SetBinContent(45,1);

      }
    }

    //cout<<"3"<<endl;

    if(R.gen_postFSR_pt->size()>=2.){   
      for(unsigned int m=0;m<R.gen_postFSR_eta->size();m++){

	if(fabs(R.gen_postFSR_eta->at(index2[m])) < 2.5 && !(fabs(R.gen_postFSR_eta->at(index2[m])) > 1.4442 && fabs(R.gen_postFSR_eta->at(index2[m])) < 1.566)){

	  gen_elec_post = gen_elec_post + 1;

	  newgPost_Pt.push_back(R.gen_postFSR_pt->at(index2[m]));
	  newgPost_Eta.push_back(R.gen_postFSR_eta->at(index2[m]));
	  newgPost_Enr.push_back(R.gen_postFSR_ene->at(index2[m]));
	  newgPost_Phi.push_back(R.gen_postFSR_phi->at(index2[m]));
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
	genPostFSR_Mass->Fill(ZMass_di_gPost,R.theWeight);
	//genPostFSR_Mass->SetBinContent(45,1);

      }
    }

    // Reco-gen matching
    if(gen_elec_post==2){
      for(unsigned int ireco = 0; ireco < neweleEta.size(); ireco++){
	double dR_comp_post = 1000.;
	for(unsigned int igen = 0; igen < newgPost_Eta.size(); igen++){
	  dR = R.deltaR(neweleEta[ireco], newelePhi[ireco], newgPost_Eta[igen], newgPost_Phi[igen]);

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

    recoPostFSRinAcc_Mass->Fill(ZMass_R_post, R.theWeight);
    genPostFSRinAcc_Mass->Fill(ZMass_G_post, R.theWeight);
    responsePost->Fill(ZMass_G_post,ZMass_R_post,R.theWeight);

    // det response
    if (insideMassRange(ZMass_R_post) || insideMassRange(ZMass_G_post)) {
      if (!insideMassRange(ZMass_R_post) &&  insideMassRange(ZMass_G_post))
      resp_DetRes->Miss(ZMass_G_post, R.theWeight);
      else if ( insideMassRange(ZMass_R_post) && !insideMassRange(ZMass_G_post))
	resp_DetRes->Fake(ZMass_R_post, R.theWeight);
      else resp_DetRes->Fill( ZMass_R_post, ZMass_G_post, R.theWeight );
    }


    // FSR response
    if (!R.tauFlag) {
      if (R.nEle>=2) {

	ele1.SetPtEtaPhiE(R.pt->at(0),R.eta->at(0),R.phi->at(0),R.energy->at(0));
	ele2.SetPtEtaPhiE(R.pt->at(1),R.eta->at(1),R.phi->at(1),R.energy->at(1));
	dielectron=ele1+ele2;
	Double_t Z_Mass_postFSR4pi = dielectron.M();
	Double_t Z_Mass_preFSR4pi = 0;

	int nPreFsr=R.gen_preFSR->GetEntriesFast();
	//std::cout << "here are nPreFsr=" << nPreFsr << " electrons\n";
	if (nPreFsr>=2) {
	  int gen_index[nPreFsr];
	  float gen_pt[nPreFsr];
	  for (int i=0; i<nPreFsr; i++) {
	    const TLorentzVector *e= (TLorentzVector*)(*R.gen_preFSR)[i];
	    std::cout << " i = " << i << ", pT=" << e->Perp() << "\n";
	    gen_pt[i]= e->Perp();
	  }
	  TMath::Sort(int(sizeof(gen_pt)/sizeof(gen_pt[0])), gen_pt,gen_index,true);

	  // Assume that 2 leading gen electrons are of our interest
	  const TLorentzVector *gen_e1= (TLorentzVector*)(*R.gen_preFSR)[gen_index[0]];
	  const TLorentzVector *gen_e2= (TLorentzVector*)(*R.gen_preFSR)[gen_index[1]];
	  std::cout << "leading electrons: " << gen_e1->Perp() << ", " << gen_e2->Perp() << "\n";
	  dielectron=(*gen_e1) + (*gen_e2);
	  Z_Mass_preFSR4pi= dielectron.M();
	}

	genPostFSR4pi_Mass->Fill(Z_Mass_postFSR4pi, R.theWeight);
	genPreFSR4pi_Mass->Fill (Z_Mass_preFSR4pi, R.theWeight);
	fsrResponse4pi->Fill(Z_Mass_preFSR4pi,Z_Mass_postFSR4pi,R.theWeight);

	// FSR response
	if (insideMassRange(Z_Mass_postFSR4pi) || insideMassRange(Z_Mass_preFSR4pi)) {
	  if (!insideMassRange(Z_Mass_postFSR4pi) &&  insideMassRange(Z_Mass_preFSR4pi))
	    resp_FSR->Miss(Z_Mass_preFSR4pi, R.theWeight);
	  else if ( insideMassRange(Z_Mass_postFSR4pi) && !insideMassRange(Z_Mass_preFSR4pi))
	    resp_FSR->Fake(Z_Mass_postFSR4pi, R.theWeight);
	  else resp_FSR->Fill( Z_Mass_postFSR4pi, Z_Mass_preFSR4pi, R.theWeight );
	}

      }
    }


  } // event

  resp_DetRes->Write();
  resp_FSR->Write();

  file->Write();
  file->Close();
  // fin.Close();
}
