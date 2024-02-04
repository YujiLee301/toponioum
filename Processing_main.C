/*
Author: Yuji Li(Fudan University)
Email: yuji.li@cern.ch
Github: YujiLee301
date: 04/02/2024
*/

#include "TLorentzVector.h"
#include "TH1F.h"
#include <vector>
#include <math.h>
#include "RooStats/RooStatsUtils.h"
#include <iostream>
using namespace RooFit;
using namespace RooStats;
using namespace std;

class EvtProducer{
    public:
      float mumuEff, bbEff;
      float muonPtcut=5;
      float jetPtCut;
      float muonEtacut=2.4; float jetEtacut=2.5;
      int Nbjetcut; //0,1,2
      float SMbbRvalue, SigbkgRvalue; // Rvalue: N(ee->bb)/N(ee->mumu) in a dataset
      float Luminosity;
      void SetInputfiles(const char* SMbbsampleName_,const char* SigbkgsampleName_,const char* mumuSampleName_,const char* bbSampleName_);
      void setTotalEvts(int nEvtbb_,int nEvtmumu_,int nEvtSig_,int nEvtSigBkg_);
      void setSM_bb_mumu_XS(float nEvtmumuXS_, float nEvtbbXS_);
      void MeasurebbEff();
      void MeasuremumuEff();
      float getRvalue(const char* SampleName);
      float GetmumuEvts(const char* inputName);
      float GetbbEvts(const char* inputName);

      float DetectedNsigbkg, DetectedNbkg;
      void getbbDetectedEvtNum();

      float DetectedNsig, JtbbExpectedXS;
      void ExtractJtbbXS();
      
      double relativeBkgUncert = 0.1;
      double pExp,zExp;

    private:
      int nEvtbb, nEvtmumu, nEvtSig, nEvtSigBkg;
      float nEvtmumuXS,nEvtbbXS;
      const char* SMbbsampleName;const char* SigbkgsampleName;const char* mumuSampleName;const char* bbSampleName;

};

void EvtProducer::SetInputfiles(const char* SMbbsampleName_,const char* SigbkgsampleName_,const char* mumuSampleName_,const char* bbSampleName_){
    SMbbsampleName = SMbbsampleName_; SigbkgsampleName = SigbkgsampleName_; 
    mumuSampleName = mumuSampleName_; bbSampleName = bbSampleName_;
}

void EvtProducer::setTotalEvts(int nEvtbb_,int nEvtmumu_,int nEvtSig_,int nEvtSigBkg_){
    nEvtbb = nEvtbb_; nEvtmumu = nEvtmumu_; nEvtSig = nEvtSig_; nEvtSigBkg = nEvtSigBkg_;
}
void EvtProducer::setSM_bb_mumu_XS(float nEvtmumuXS_, float nEvtbbXS_){
    nEvtmumuXS = nEvtmumuXS_; nEvtbbXS = nEvtbbXS_;
}
void EvtProducer::MeasurebbEff(){
    float Nevts;
    Nevts = GetbbEvts(bbSampleName);
    bbEff = Nevts/nEvtbb;
}
void EvtProducer::MeasuremumuEff(){
    float Nevts;
    Nevts = GetmumuEvts(mumuSampleName);
    mumuEff = Nevts/nEvtmumu;
}
float EvtProducer::getRvalue(const char* SampleName){
    float Nevtsbb, NevtsMumu;
    Nevtsbb = GetbbEvts(SampleName);
    NevtsMumu = GetmumuEvts(SampleName);
    return Nevtsbb/NevtsMumu;
}
float EvtProducer::GetmumuEvts(const char* inputName){
    float Nevts=0;
    TFile * inFile = TFile :: Open(inputName, "READ");
    TTree * tree = (TTree*) inFile->Get("Delphes");
    int Muon_size;
  
    Float_t* Muon_pt = new Float_t[10];
    Float_t* Muon_eta = new Float_t[10];
    Float_t* Muon_phi = new Float_t[10];

    tree->SetBranchAddress("Muon.PT", Muon_pt);
    tree->SetBranchAddress("Muon.Eta", Muon_eta);
    tree->SetBranchAddress("Muon.Phi", Muon_phi);
    tree->SetBranchAddress("Muon_size", &Muon_size);

    Long64_t totalEntry = tree->GetEntries();
  
    for ( Long64_t entry = 0; entry < totalEntry; ++entry)
    { 
        tree -> GetEntry(entry);
        if (Muon_size < 2) continue;
        int Nmuon = 0;
        for(int k=0; k<Muon_size;k++){
            if((Muon_pt[k]>muonPtcut)&&(abs(Muon_eta[k])<muonEtacut)){
            Nmuon++;
            } 
        } 
        if(Nmuon<2) continue;
        Nevts++;
    }
    inFile->Close();
    return Nevts;
}
float EvtProducer::GetbbEvts(const char* inputName){
    float Nevts=0;
    TFile * inFile = TFile :: Open(inputName, "READ");
    TTree * tree = (TTree*) inFile->Get("Delphes");
    int Jet_size;
  
    Float_t* Jet_pt = new Float_t[10];
    Float_t* Jet_eta = new Float_t[10];
    Float_t* Jet_phi = new Float_t[10];
    Float_t* Jet_m = new Float_t[10];
    UInt_t* Jet_BTag = new UInt_t[10];

    tree->SetBranchAddress("Jet.PT", Jet_pt);
    tree->SetBranchAddress("Jet.Eta", Jet_eta);
    tree->SetBranchAddress("Jet.Phi", Jet_phi);
    tree->SetBranchAddress("Jet.Mass", Jet_m);
    tree->SetBranchAddress("Jet_size", &Jet_size);
    tree->SetBranchAddress("Jet.BTag", Jet_BTag);

    Long64_t totalEntry = tree->GetEntries();
  
    for ( Long64_t entry = 0; entry < totalEntry; ++entry)
    { 
        tree -> GetEntry(entry);
        if (Jet_size < 2) continue;
        int Njet = 0;
        int bjetcount = 0;
        vector<int> jetidx,jetidx_;
        for(int k=0; k<Jet_size;k++){
            if(Jet_BTag[k]) bjetcount++;
            if(Jet_pt[k]>jetPtCut){
            Njet++;
            } 
        } 
        if(Njet<2) continue;
        if(bjetcount<Nbjetcut) continue;
        Nevts++;
    }
    bbEff = Nevts/nEvtbb;
    inFile->Close();
    return Nevts;
}

void EvtProducer::getbbDetectedEvtNum(){
    float Nmumu;
    Nmumu = Luminosity*nEvtmumuXS*mumuEff;
    DetectedNsigbkg = Nmumu * getRvalue(SigbkgsampleName);
    DetectedNbkg = Nmumu * getRvalue(SMbbsampleName);
}

void EvtProducer::ExtractJtbbXS(){
    DetectedNsig = DetectedNsigbkg - DetectedNbkg;
    JtbbExpectedXS = DetectedNsig/(Luminosity*bbEff);
    relativeBkgUncert = 1/sqrt(DetectedNbkg);
    pExp = NumberCountingUtils::BinomialExpP(DetectedNsig, DetectedNbkg, relativeBkgUncert);
    zExp = NumberCountingUtils::BinomialExpZ(DetectedNsig, DetectedNbkg, relativeBkgUncert);
}

class InputTags{
    public:
      EvtProducer Ana;
      float Energy = 345.91;
      float jetPtCut = 150;
      float Nbjetcut = 1;
      float Luminosity = 200; 
      float totalGENevts = 100000;
      float bbXS = 826.21;
      float mumuXS = 954.11;
      const char* SMbbsampleName;const char* SigbkgsampleName;const char* mumuSampleName;const char* bbSampleName;
};
void Processing(InputTags Input){
    Input.Ana.SetInputfiles("delphes_output_eebb_B_345p91.root","delphes_output_eebb_SB_345p91.root","delphes_output_eemumu_345p91.root","delphes_output_eebb_345p91.root");
    Input.Ana.jetPtCut = Input.jetPtCut; Input.Ana.Nbjetcut = Input.Nbjetcut; 
    Input.Ana.Luminosity = Input.Luminosity;
    Input.Ana.setTotalEvts(Input.totalGENevts,Input.totalGENevts,Input.totalGENevts,Input.totalGENevts);
    Input.Ana.setSM_bb_mumu_XS(Input.mumuXS, Input.bbXS); //[fb]
    Input.Ana.MeasurebbEff();
    Input.Ana.MeasuremumuEff();
    Input.Ana.getbbDetectedEvtNum();
    Input.Ana.ExtractJtbbXS();
    cout<<"SigNum:"<<Input.Ana.DetectedNsig<<" S+BNum:"<<Input.Ana.DetectedNsigbkg<< "  Gaussian sigma = "<< -1*Input.Ana.zExp << endl;
}


void Processing_main(){
    InputTags input;
    input.Energy = 345.91;
    input.jetPtCut = 150;
    input.Nbjetcut = 1;
    input.Luminosity = 200; 
    input.totalGENevts = 100000;
    input.bbXS = 826.21;//[fb]
    input.mumuXS = 954.11;//[fb]
    input.SMbbsampleName = "delphes_output_eebb_B_345p91.root";
    input.SigbkgsampleName = "delphes_output_eebb_SB_345p91.root";
    input.mumuSampleName = "delphes_output_eemumu_345p91.root";
    input.bbSampleName = "delphes_output_eebb_345p91.root";
    Processing(input);
}

/*
All the input files could be downloaded from this link:
https://cernbox.cern.ch/

Out-of-box Example:
root -q Processing_main.C
*/