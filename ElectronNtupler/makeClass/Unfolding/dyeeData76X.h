//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 21 12:48:51 2016 by ROOT version 6.02/05
// from TTree ElectronTree/Electron data
// found on file: /tmp/andriusj/DY_50toInf.root
//////////////////////////////////////////////////////////

#ifndef dyeeData76X_h
#define dyeeData76X_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <sstream>


class dyeeData76X {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        RunNo;
   Double_t        EvtNo;
   Double_t        Lumi;
   Double_t        Bunch;
   Double_t        theWeight;
   Int_t           pvNTracks;
   Int_t           nPV;
   Int_t           nPU;
   Int_t           nPUTrue;
   Float_t         rho;
   Bool_t          Ele23_WPLoose;
   Bool_t          Ele27_WP85;
   Bool_t          IsoMu20;
   Bool_t          Ele17_Ele12;
   Bool_t          Mu8_Ele17;
   Bool_t          Mu17_Ele12;
   Bool_t          Mu8_Ele23;
   vector<double>  *etPhoton;
   Int_t           singlePhoton;
   Int_t           prescalePhoton;
   vector<double>  *pt_Ele23;
   vector<double>  *eta_Ele23;
   vector<double>  *phi_Ele23;
   Int_t           nEle;
   Int_t           nGenEle;
   Int_t           nGenTau;
   Int_t           tauFlag;
   vector<bool>    *fromHProcessFinalState;
   vector<bool>    *fromHProcessDecayed;
   vector<float>   *genPhoton_Px;
   vector<float>   *genPhoton_Py;
   vector<float>   *genPhoton_Pz;
   vector<float>   *genPhoton_Pt;
   vector<float>   *genPhoton_Eta;
   vector<float>   *genPhoton_Rap;
   vector<float>   *genPhoton_Phi;
   vector<float>   *genPhoton_En;
   vector<float>   *genPostFSR_Px;
   vector<float>   *genPostFSR_Py;
   vector<float>   *genPostFSR_Pz;
   vector<float>   *genPostFSR_Pt;
   vector<float>   *genPostFSR_Eta;
   vector<float>   *genPostFSR_Rap;
   vector<float>   *genPostFSR_Phi;
   vector<float>   *genPostFSR_En;
   vector<float>   *genPreFSR_Px;
   vector<float>   *genPreFSR_Py;
   vector<float>   *genPreFSR_Pz;
   vector<float>   *genPreFSR_Pt;
   vector<float>   *genPreFSR_Eta;
   vector<float>   *genPreFSR_Rap;
   vector<float>   *genPreFSR_Phi;
   vector<float>   *genPreFSR_En;
   vector<float>   *ptElec;
   vector<float>   *etaElec;
   vector<float>   *rapElec;
   vector<float>   *phiElec;
   vector<float>   *energyElec;
   vector<float>   *massElec;
   vector<float>   *chargeElec;
   vector<float>   *enSC;
   vector<float>   *preEnSC;
   vector<float>   *rawEnSC;
   vector<float>   *etSC;
   vector<float>   *etaSC;
   vector<float>   *phiSC;
   vector<float>   *full5x5_sigmaIetaIeta;
   vector<float>   *E1x5;
   vector<float>   *E2x5;
   vector<float>   *E5x5;
   vector<float>   *hOverE;
   vector<float>   *etaScWidth;
   vector<float>   *phiScWidth;
   vector<float>   *r9;
   vector<float>   *dEtaIn;
   vector<float>   *dPhiIn;
   vector<float>   *isoChargedHadrons;
   vector<float>   *isoNeutralHadrons;
   vector<float>   *isoPhotons;
   vector<float>   *isoChargedFromPU;
   vector<float>   *isoDeltaBeta;
   vector<float>   *isoRho;
   vector<float>   *ooEmooP;
   vector<float>   *d0;
   vector<float>   *dz;
   vector<int>     *expectedMissingInnerHits;
   vector<int>     *passConversionVeto;
   vector<float>   *brem;
   vector<int>     *isTrue;
   vector<int>     *passVetoId;
   vector<int>     *passLooseId;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
   vector<int>     *passHEEPId;
   vector<int>     *isPassMedium_NoPt;
   vector<int>     *isPassMedium_NoScEta;
   vector<int>     *isPassMedium_NoDEta;
   vector<int>     *isPassMedium_NoDPhi;
   vector<int>     *isPassMedium_NoSigmaEtaEta;
   vector<int>     *isPassMedium_NoHOverE;
   vector<int>     *isPassMedium_NoDxy;
   vector<int>     *isPassMedium_NoDz;
   vector<int>     *isPassMedium_NoEInvP;
   vector<int>     *isPassMedium_NoPFIso;
   vector<int>     *isPassMedium_NoConVeto;
   vector<int>     *isPassMedium_NoMissHits;
   vector<int>     *eleEcalDrivenSeed;
   vector<float>   *eleInBarrel;
   vector<float>   *eleInEndcap;
   Int_t           nMuons;
   vector<bool>    *isLooseMuon;
   vector<bool>    *isTightMuon;
   vector<bool>    *isHighPtMuon;
   vector<double>  *ptMuon;
   vector<double>  *etaMuon;
   vector<double>  *phiMuon;
   vector<double>  *energyMuon;
   vector<double>  *chargeMuon;
   vector<double>  *isoChargedHadronPfR04Muon;
   vector<double>  *isoNeutralHadronPfR04Muon;
   vector<double>  *isoGammaPfR04Muon;
   vector<double>  *isoChargedFromPUMuon;
   vector<double>  *isoPFMuon;
   vector<double>  *isoTrkMuon;
   vector<double>  *metPt;
   vector<double>  *metPhi;
   vector<double>  *metSumEt;
   Int_t           nPhotons;
   vector<float>   *ptPhoton;
   vector<float>   *etaPhoton;
   vector<float>   *phiPhoton;

   // List of branches
   TBranch        *b_RunNo;   //!
   TBranch        *b_EvtNo;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_Bunch;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_pvNTracks;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nPUTrue;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_Ele23_WPLoose;   //!
   TBranch        *b_Ele27_WP85;   //!
   TBranch        *b_IsoMu20;   //!
   TBranch        *b_Ele17_Ele12;   //!
   TBranch        *b_Mu8_Ele17;   //!
   TBranch        *b_Mu17_Ele12;   //!
   TBranch        *b_Mu8_Ele23;   //!
   TBranch        *b_etPhoton;   //!
   TBranch        *b_singlePhoton;   //!
   TBranch        *b_prescalePhoton;   //!
   TBranch        *b_pt_Ele23;   //!
   TBranch        *b_eta_Ele23;   //!
   TBranch        *b_phi_Ele23;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_nGenEle;   //!
   TBranch        *b_nGenTau;   //!
   TBranch        *b_tauFlag;   //!
   TBranch        *b_fromHProcessFinalState;   //!
   TBranch        *b_fromHProcessDecayed;   //!
   TBranch        *b_genPhoton_Px;   //!
   TBranch        *b_genPhoton_Py;   //!
   TBranch        *b_genPhoton_Pz;   //!
   TBranch        *b_genPhoton_Pt;   //!
   TBranch        *b_genPhoton_Eta;   //!
   TBranch        *b_genPhoton_Rap;   //!
   TBranch        *b_genPhoton_Phi;   //!
   TBranch        *b_genPhoton_En;   //!
   TBranch        *b_genPostFSR_Px;   //!
   TBranch        *b_genPostFSR_Py;   //!
   TBranch        *b_genPostFSR_Pz;   //!
   TBranch        *b_genPostFSR_Pt;   //!
   TBranch        *b_genPostFSR_Eta;   //!
   TBranch        *b_genPostFSR_Rap;   //!
   TBranch        *b_genPostFSR_Phi;   //!
   TBranch        *b_genPostFSR_En;   //!
   TBranch        *b_genPreFSR_Px;   //!
   TBranch        *b_genPreFSR_Py;   //!
   TBranch        *b_genPreFSR_Pz;   //!
   TBranch        *b_genPreFSR_Pt;   //!
   TBranch        *b_genPreFSR_Eta;   //!
   TBranch        *b_genPreFSR_Rap;   //!
   TBranch        *b_genPreFSR_Phi;   //!
   TBranch        *b_genPreFSR_En;   //!
   TBranch        *b_ptElec;   //!
   TBranch        *b_etaElec;   //!
   TBranch        *b_rapElec;   //!
   TBranch        *b_phiElec;   //!
   TBranch        *b_energyElec;   //!
   TBranch        *b_massElec;   //!
   TBranch        *b_chargeElec;   //!
   TBranch        *b_enSC;   //!
   TBranch        *b_preEnSC;   //!
   TBranch        *b_rawEnSC;   //!
   TBranch        *b_etSC;   //!
   TBranch        *b_etaSC;   //!
   TBranch        *b_phiSC;   //!
   TBranch        *b_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_E1x5;   //!
   TBranch        *b_E2x5;   //!
   TBranch        *b_E5x5;   //!
   TBranch        *b_hOverE;   //!
   TBranch        *b_etaScWidth;   //!
   TBranch        *b_phiScWidth;   //!
   TBranch        *b_r9;   //!
   TBranch        *b_dEtaIn;   //!
   TBranch        *b_dPhiIn;   //!
   TBranch        *b_isoChargedHadrons;   //!
   TBranch        *b_isoNeutralHadrons;   //!
   TBranch        *b_isoPhotons;   //!
   TBranch        *b_isoChargedFromPU;   //!
   TBranch        *b_isoDeltaBeta;   //!
   TBranch        *b_isoRho;   //!
   TBranch        *b_ooEmooP;   //!
   TBranch        *b_d0;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_expectedMissingInnerHits;   //!
   TBranch        *b_passConversionVeto;   //!
   TBranch        *b_brem;   //!
   TBranch        *b_isTrue;   //!
   TBranch        *b_passVetoId;   //!
   TBranch        *b_passLooseId;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_passTightId;   //!
   TBranch        *b_passHEEPId;   //!
   TBranch        *b_isPassMedium_NoPt;   //!
   TBranch        *b_isPassMedium_NoScEta;   //!
   TBranch        *b_isPassMedium_NoDEta;   //!
   TBranch        *b_isPassMedium_NoDPhi;   //!
   TBranch        *b_isPassMedium_NoSigmaEtaEta;   //!
   TBranch        *b_isPassMedium_NoHOverE;   //!
   TBranch        *b_isPassMedium_NoDxy;   //!
   TBranch        *b_isPassMedium_NoDz;   //!
   TBranch        *b_isPassMedium_NoEInvP;   //!
   TBranch        *b_isPassMedium_NoPFIso;   //!
   TBranch        *b_isPassMedium_NoConVeto;   //!
   TBranch        *b_isPassMedium_NoMissHits;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_eleInBarrel;   //!
   TBranch        *b_eleInEndcap;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_isLooseMuon;   //!
   TBranch        *b_isTightMuon;   //!
   TBranch        *b_isHighPtMuon;   //!
   TBranch        *b_ptMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_isoChargedHadronPfR04Muon;   //!
   TBranch        *b_isoNeutralHadronPfR04Muon;   //!
   TBranch        *b_isoGammaPfR04Muon;   //!
   TBranch        *b_isoChargedFromPUMuon;   //!
   TBranch        *b_isoPFMuon;   //!
   TBranch        *b_isoTrkMuon;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_metSumEt;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_ptPhoton;   //!
   TBranch        *b_etaPhoton;   //!
   TBranch        *b_phiPhoton;   //!

   dyeeData76X(TString inpFiles="", TString inpTree="ntupler/ElectronTree");
   virtual ~dyeeData76X();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual int      Init(TString inpFiles);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   void DeactivateBranches() { fChain->SetBranchStatus("*",0); }

   void ActivateBranches(TString brList) {
     TString brName;
     std::stringstream ss(brList.Data());
     while (!ss.eof()) {
       ss >> brName;
       if (brName.Length()) {
	 std::cout << "activate branch " << brName << "\n";
	 fChain->SetBranchStatus(brName,1);
       }
     }
   }

   inline static
   Double_t deltaPhi(Double_t phi1, Double_t phi2)
     {
       Double_t pi = 3.1415927;
       Double_t dphi = fabs(phi1 - phi2);
       if(dphi >= pi) dphi = 2. * pi - dphi;
       return dphi;
     }

   inline static
   Double_t deltaEta(Double_t eta1, Double_t eta2)
     {
       Double_t deta = fabs(eta1 - eta2);
       return deta;
     }

   inline static
   Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
     {
       Double_t deta = deltaEta(eta1, eta2);
       Double_t dphi = deltaPhi(phi1, phi2);
       Double_t dr = sqrt(deta*deta + dphi*dphi);
       return dr;
     }


};

#endif

#ifdef dyeeData76X_cxx
dyeeData76X::dyeeData76X(TString inpFiles, TString treeName) :
  fChain(new TChain(treeName))
{
  if (!Init(inpFiles)) std::cout << "dyeeData76X creation error\n";
}

dyeeData76X::~dyeeData76X()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dyeeData76X::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dyeeData76X::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

int dyeeData76X::Init(TString inpFiles)
{

  if (inpFiles.Length()==0) return 0;
  TString fn;
  std::stringstream ss(inpFiles.Data());
  while (!ss.eof()) {
    ss >> fn;
    if (fn.Length()) {
      std::cout << "adding file <" << fn << ">\n";
      fChain->Add(fn);
    }
  }

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   etPhoton = 0;
   pt_Ele23 = 0;
   eta_Ele23 = 0;
   phi_Ele23 = 0;
   fromHProcessFinalState = 0;
   fromHProcessDecayed = 0;
   genPhoton_Px = 0;
   genPhoton_Py = 0;
   genPhoton_Pz = 0;
   genPhoton_Pt = 0;
   genPhoton_Eta = 0;
   genPhoton_Rap = 0;
   genPhoton_Phi = 0;
   genPhoton_En = 0;
   genPostFSR_Px = 0;
   genPostFSR_Py = 0;
   genPostFSR_Pz = 0;
   genPostFSR_Pt = 0;
   genPostFSR_Eta = 0;
   genPostFSR_Rap = 0;
   genPostFSR_Phi = 0;
   genPostFSR_En = 0;
   genPreFSR_Px = 0;
   genPreFSR_Py = 0;
   genPreFSR_Pz = 0;
   genPreFSR_Pt = 0;
   genPreFSR_Eta = 0;
   genPreFSR_Rap = 0;
   genPreFSR_Phi = 0;
   genPreFSR_En = 0;
   ptElec = 0;
   etaElec = 0;
   rapElec = 0;
   phiElec = 0;
   energyElec = 0;
   massElec = 0;
   chargeElec = 0;
   enSC = 0;
   preEnSC = 0;
   rawEnSC = 0;
   etSC = 0;
   etaSC = 0;
   phiSC = 0;
   full5x5_sigmaIetaIeta = 0;
   E1x5 = 0;
   E2x5 = 0;
   E5x5 = 0;
   hOverE = 0;
   etaScWidth = 0;
   phiScWidth = 0;
   r9 = 0;
   dEtaIn = 0;
   dPhiIn = 0;
   isoChargedHadrons = 0;
   isoNeutralHadrons = 0;
   isoPhotons = 0;
   isoChargedFromPU = 0;
   isoDeltaBeta = 0;
   isoRho = 0;
   ooEmooP = 0;
   d0 = 0;
   dz = 0;
   expectedMissingInnerHits = 0;
   passConversionVeto = 0;
   brem = 0;
   isTrue = 0;
   passVetoId = 0;
   passLooseId = 0;
   passMediumId = 0;
   passTightId = 0;
   passHEEPId = 0;
   isPassMedium_NoPt = 0;
   isPassMedium_NoScEta = 0;
   isPassMedium_NoDEta = 0;
   isPassMedium_NoDPhi = 0;
   isPassMedium_NoSigmaEtaEta = 0;
   isPassMedium_NoHOverE = 0;
   isPassMedium_NoDxy = 0;
   isPassMedium_NoDz = 0;
   isPassMedium_NoEInvP = 0;
   isPassMedium_NoPFIso = 0;
   isPassMedium_NoConVeto = 0;
   isPassMedium_NoMissHits = 0;
   eleEcalDrivenSeed = 0;
   eleInBarrel = 0;
   eleInEndcap = 0;
   isLooseMuon = 0;
   isTightMuon = 0;
   isHighPtMuon = 0;
   ptMuon = 0;
   etaMuon = 0;
   phiMuon = 0;
   energyMuon = 0;
   chargeMuon = 0;
   isoChargedHadronPfR04Muon = 0;
   isoNeutralHadronPfR04Muon = 0;
   isoGammaPfR04Muon = 0;
   isoChargedFromPUMuon = 0;
   isoPFMuon = 0;
   isoTrkMuon = 0;
   metPt = 0;
   metPhi = 0;
   metSumEt = 0;
   ptPhoton = 0;
   etaPhoton = 0;
   phiPhoton = 0;
   // Set branch addresses and branch pointers
   //if (!tree) return;
   //fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("EvtNo", &EvtNo, &b_EvtNo);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("Bunch", &Bunch, &b_Bunch);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("pvNTracks", &pvNTracks, &b_pvNTracks);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nPUTrue", &nPUTrue, &b_nPUTrue);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose, &b_Ele23_WPLoose);
   fChain->SetBranchAddress("Ele27_WP85", &Ele27_WP85, &b_Ele27_WP85);
   fChain->SetBranchAddress("IsoMu20", &IsoMu20, &b_IsoMu20);
   fChain->SetBranchAddress("Ele17_Ele12", &Ele17_Ele12, &b_Ele17_Ele12);
   fChain->SetBranchAddress("Mu8_Ele17", &Mu8_Ele17, &b_Mu8_Ele17);
   fChain->SetBranchAddress("Mu17_Ele12", &Mu17_Ele12, &b_Mu17_Ele12);
   fChain->SetBranchAddress("Mu8_Ele23", &Mu8_Ele23, &b_Mu8_Ele23);
   fChain->SetBranchAddress("etPhoton", &etPhoton, &b_etPhoton);
   fChain->SetBranchAddress("singlePhoton", &singlePhoton, &b_singlePhoton);
   fChain->SetBranchAddress("prescalePhoton", &prescalePhoton, &b_prescalePhoton);
   fChain->SetBranchAddress("pt_Ele23", &pt_Ele23, &b_pt_Ele23);
   fChain->SetBranchAddress("eta_Ele23", &eta_Ele23, &b_eta_Ele23);
   fChain->SetBranchAddress("phi_Ele23", &phi_Ele23, &b_phi_Ele23);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("nGenEle", &nGenEle, &b_nGenEle);
   fChain->SetBranchAddress("nGenTau", &nGenTau, &b_nGenTau);
   fChain->SetBranchAddress("tauFlag", &tauFlag, &b_tauFlag);
   fChain->SetBranchAddress("fromHProcessFinalState", &fromHProcessFinalState, &b_fromHProcessFinalState);
   fChain->SetBranchAddress("fromHProcessDecayed", &fromHProcessDecayed, &b_fromHProcessDecayed);
   fChain->SetBranchAddress("genPhoton_Px", &genPhoton_Px, &b_genPhoton_Px);
   fChain->SetBranchAddress("genPhoton_Py", &genPhoton_Py, &b_genPhoton_Py);
   fChain->SetBranchAddress("genPhoton_Pz", &genPhoton_Pz, &b_genPhoton_Pz);
   fChain->SetBranchAddress("genPhoton_Pt", &genPhoton_Pt, &b_genPhoton_Pt);
   fChain->SetBranchAddress("genPhoton_Eta", &genPhoton_Eta, &b_genPhoton_Eta);
   fChain->SetBranchAddress("genPhoton_Rap", &genPhoton_Rap, &b_genPhoton_Rap);
   fChain->SetBranchAddress("genPhoton_Phi", &genPhoton_Phi, &b_genPhoton_Phi);
   fChain->SetBranchAddress("genPhoton_En", &genPhoton_En, &b_genPhoton_En);
   fChain->SetBranchAddress("genPostFSR_Px", &genPostFSR_Px, &b_genPostFSR_Px);
   fChain->SetBranchAddress("genPostFSR_Py", &genPostFSR_Py, &b_genPostFSR_Py);
   fChain->SetBranchAddress("genPostFSR_Pz", &genPostFSR_Pz, &b_genPostFSR_Pz);
   fChain->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt, &b_genPostFSR_Pt);
   fChain->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta, &b_genPostFSR_Eta);
   fChain->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap, &b_genPostFSR_Rap);
   fChain->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi, &b_genPostFSR_Phi);
   fChain->SetBranchAddress("genPostFSR_En", &genPostFSR_En, &b_genPostFSR_En);
   fChain->SetBranchAddress("genPreFSR_Px", &genPreFSR_Px, &b_genPreFSR_Px);
   fChain->SetBranchAddress("genPreFSR_Py", &genPreFSR_Py, &b_genPreFSR_Py);
   fChain->SetBranchAddress("genPreFSR_Pz", &genPreFSR_Pz, &b_genPreFSR_Pz);
   fChain->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt, &b_genPreFSR_Pt);
   fChain->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta, &b_genPreFSR_Eta);
   fChain->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap, &b_genPreFSR_Rap);
   fChain->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi, &b_genPreFSR_Phi);
   fChain->SetBranchAddress("genPreFSR_En", &genPreFSR_En, &b_genPreFSR_En);
   fChain->SetBranchAddress("ptElec", &ptElec, &b_ptElec);
   fChain->SetBranchAddress("etaElec", &etaElec, &b_etaElec);
   fChain->SetBranchAddress("rapElec", &rapElec, &b_rapElec);
   fChain->SetBranchAddress("phiElec", &phiElec, &b_phiElec);
   fChain->SetBranchAddress("energyElec", &energyElec, &b_energyElec);
   fChain->SetBranchAddress("massElec", &massElec, &b_massElec);
   fChain->SetBranchAddress("chargeElec", &chargeElec, &b_chargeElec);
   fChain->SetBranchAddress("enSC", &enSC, &b_enSC);
   fChain->SetBranchAddress("preEnSC", &preEnSC, &b_preEnSC);
   fChain->SetBranchAddress("rawEnSC", &rawEnSC, &b_rawEnSC);
   fChain->SetBranchAddress("etSC", &etSC, &b_etSC);
   fChain->SetBranchAddress("etaSC", &etaSC, &b_etaSC);
   fChain->SetBranchAddress("phiSC", &phiSC, &b_phiSC);
   fChain->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta, &b_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("E1x5", &E1x5, &b_E1x5);
   fChain->SetBranchAddress("E2x5", &E2x5, &b_E2x5);
   fChain->SetBranchAddress("E5x5", &E5x5, &b_E5x5);
   fChain->SetBranchAddress("hOverE", &hOverE, &b_hOverE);
   fChain->SetBranchAddress("etaScWidth", &etaScWidth, &b_etaScWidth);
   fChain->SetBranchAddress("phiScWidth", &phiScWidth, &b_phiScWidth);
   fChain->SetBranchAddress("r9", &r9, &b_r9);
   fChain->SetBranchAddress("dEtaIn", &dEtaIn, &b_dEtaIn);
   fChain->SetBranchAddress("dPhiIn", &dPhiIn, &b_dPhiIn);
   fChain->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
   fChain->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
   fChain->SetBranchAddress("isoPhotons", &isoPhotons, &b_isoPhotons);
   fChain->SetBranchAddress("isoChargedFromPU", &isoChargedFromPU, &b_isoChargedFromPU);
   fChain->SetBranchAddress("isoDeltaBeta", &isoDeltaBeta, &b_isoDeltaBeta);
   fChain->SetBranchAddress("isoRho", &isoRho, &b_isoRho);
   fChain->SetBranchAddress("ooEmooP", &ooEmooP, &b_ooEmooP);
   fChain->SetBranchAddress("d0", &d0, &b_d0);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("expectedMissingInnerHits", &expectedMissingInnerHits, &b_expectedMissingInnerHits);
   fChain->SetBranchAddress("passConversionVeto", &passConversionVeto, &b_passConversionVeto);
   fChain->SetBranchAddress("brem", &brem, &b_brem);
   fChain->SetBranchAddress("isTrue", &isTrue, &b_isTrue);
   fChain->SetBranchAddress("passVetoId", &passVetoId, &b_passVetoId);
   fChain->SetBranchAddress("passLooseId", &passLooseId, &b_passLooseId);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("passTightId", &passTightId, &b_passTightId);
   fChain->SetBranchAddress("passHEEPId", &passHEEPId, &b_passHEEPId);
   fChain->SetBranchAddress("isPassMedium_NoPt", &isPassMedium_NoPt, &b_isPassMedium_NoPt);
   fChain->SetBranchAddress("isPassMedium_NoScEta", &isPassMedium_NoScEta, &b_isPassMedium_NoScEta);
   fChain->SetBranchAddress("isPassMedium_NoDEta", &isPassMedium_NoDEta, &b_isPassMedium_NoDEta);
   fChain->SetBranchAddress("isPassMedium_NoDPhi", &isPassMedium_NoDPhi, &b_isPassMedium_NoDPhi);
   fChain->SetBranchAddress("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta, &b_isPassMedium_NoSigmaEtaEta);
   fChain->SetBranchAddress("isPassMedium_NoHOverE", &isPassMedium_NoHOverE, &b_isPassMedium_NoHOverE);
   fChain->SetBranchAddress("isPassMedium_NoDxy", &isPassMedium_NoDxy, &b_isPassMedium_NoDxy);
   fChain->SetBranchAddress("isPassMedium_NoDz", &isPassMedium_NoDz, &b_isPassMedium_NoDz);
   fChain->SetBranchAddress("isPassMedium_NoEInvP", &isPassMedium_NoEInvP, &b_isPassMedium_NoEInvP);
   fChain->SetBranchAddress("isPassMedium_NoPFIso", &isPassMedium_NoPFIso, &b_isPassMedium_NoPFIso);
   fChain->SetBranchAddress("isPassMedium_NoConVeto", &isPassMedium_NoConVeto, &b_isPassMedium_NoConVeto);
   fChain->SetBranchAddress("isPassMedium_NoMissHits", &isPassMedium_NoMissHits, &b_isPassMedium_NoMissHits);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleInBarrel", &eleInBarrel, &b_eleInBarrel);
   fChain->SetBranchAddress("eleInEndcap", &eleInEndcap, &b_eleInEndcap);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("isLooseMuon", &isLooseMuon, &b_isLooseMuon);
   fChain->SetBranchAddress("isTightMuon", &isTightMuon, &b_isTightMuon);
   fChain->SetBranchAddress("isHighPtMuon", &isHighPtMuon, &b_isHighPtMuon);
   fChain->SetBranchAddress("ptMuon", &ptMuon, &b_ptMuon);
   fChain->SetBranchAddress("etaMuon", &etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", &phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("energyMuon", &energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("chargeMuon", &chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("isoChargedHadronPfR04Muon", &isoChargedHadronPfR04Muon, &b_isoChargedHadronPfR04Muon);
   fChain->SetBranchAddress("isoNeutralHadronPfR04Muon", &isoNeutralHadronPfR04Muon, &b_isoNeutralHadronPfR04Muon);
   fChain->SetBranchAddress("isoGammaPfR04Muon", &isoGammaPfR04Muon, &b_isoGammaPfR04Muon);
   fChain->SetBranchAddress("isoChargedFromPUMuon", &isoChargedFromPUMuon, &b_isoChargedFromPUMuon);
   fChain->SetBranchAddress("isoPFMuon", &isoPFMuon, &b_isoPFMuon);
   fChain->SetBranchAddress("isoTrkMuon", &isoTrkMuon, &b_isoTrkMuon);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("ptPhoton", &ptPhoton, &b_ptPhoton);
   fChain->SetBranchAddress("etaPhoton", &etaPhoton, &b_etaPhoton);
   fChain->SetBranchAddress("phiPhoton", &phiPhoton, &b_phiPhoton);
   Notify();
   return 1;
}

Bool_t dyeeData76X::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dyeeData76X::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dyeeData76X::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#undef dyeeData76X_cxx
#endif // #ifdef dyeeData76X_cxx
