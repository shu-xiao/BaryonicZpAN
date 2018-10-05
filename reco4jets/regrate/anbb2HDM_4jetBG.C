#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH2F.h>
#include <TFile.h>
#include <TPaveStats.h>
#include "../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include <TLegend.h>
#include "string"
#include "../setNCUStyle.C"

#define saveEffi        true
#define doGenMatch      0
#define saveTree        1
#define corr            0

# if doGenMatch != 0 || corr != 0 
#include "../../gen2HDMsample/genMatch.C"
#endif
bool isOverlap(int i1, int i2, int j1, int j2) {
    if (i1==j1||i2==j2) return true;
    else if (i1==j2||i2==j1) return true;
    else return false;
}
void arrangeIndex(int& ind0, int& ind1, int& ind2, int& ind3) {
    bool b = false;
    if (ind0==ind1) return;
    for (int k=0;k<4;k++) {
        if (k==ind0||k==ind1) continue;
        if (!b) {
            ind2 = k;
            b = true;
        }
        else ind3 = k;
    }

}
int searchOrder(int n,int *list) {
    int k = 0;
    for (int i=0;i<4;i++) {
        if (list[n]>list[i]) k++;
    }
    return k;
}
bool isCorrect(TLorentzVector* fourJet[],TLorentzVector* genPar[]) {
    bool cor[4] = {0,0,0,0};
    const float conSize = 0.4;
    if (fourJet[0]->DeltaR(*genPar[0])<conSize&&fourJet[1]->DeltaR(*genPar[1])<conSize) cor[0] = true;
    if (fourJet[0]->DeltaR(*genPar[1])<conSize&&fourJet[1]->DeltaR(*genPar[0])<conSize) cor[1] = true;
    if (fourJet[2]->DeltaR(*genPar[2])<conSize&&fourJet[3]->DeltaR(*genPar[3])<conSize) cor[2] = true;
    if (fourJet[2]->DeltaR(*genPar[3])<conSize&&fourJet[3]->DeltaR(*genPar[2])<conSize) cor[3] = true;
    
    return (cor[0]||cor[1])&&(cor[2]||cor[3]);
}
void efferr(float nsig,float ntotal,float factor=1)
{
    float eff = nsig/ntotal;
    float err = sqrt( (1-eff)*eff/ntotal);
    cout << "efficiency = " << eff*factor << " +- " << err*factor << endl;
    //ofstream myfile;
    //myfile.open();
    //myfile << eff*factor << endl;
    //myfile << err*factor;
}
inline float caldePhi(float phi1, float phi2) {
    float dePhi = 0;
    if (abs(phi1-phi2)>TMath::Pi()) dePhi = 2*TMath::Pi() - abs(phi1-phi2);
    else dePhi = abs(phi2-phi2);
    return dePhi;
}
inline float ptAssymetry(TLorentzVector* j1, TLorentzVector* j2) {
    float minPt = min(j1->Pt(),j2->Pt());
    float deR = j1->DeltaR(*j2);
    float mj = (*j1+*j2).M();
    return pow(minPt*deR/mj,2);
}
inline float softDropAs(TLorentzVector *j1, TLorentzVector *j2, float r0 = 0.4, float beta = 0) {
    float minPt = (j1->Pt()>j2->Pt())? j2->Pt():j1->Pt();
    return minPt/(j1->Pt()+j2->Pt())*pow(r0/j1->DeltaR(*j2),beta);
}
void savenPass(int nPass[],string fileName) {
    fstream fp;
    fp.open(fileName.data(), ios::out);
    for (int i=0;i<20;i++) {
        if (nPass[i]==0) break;
        fp << (float)nPass[i] << endl;
    }
    fp.close();

}
void savenPass(int nPass[],int nPass2[],string fileName) {
    fstream fp;
    fp.open(fileName.data(), ios::out);
    for (int i=0;i<20;i++) {
        //cout << "tt" << "\t" << nPass[i] << "\t" << nPass2[i] << endl;
        if (nPass[i]==0) break;
        fp << (float)nPass[i] <<"\t"<< (float)nPass2[i] << endl;
    }
    fp.close();

}
using namespace std;
void anbb2HDM_4jetBG(int w=0, std::string inputFile="2HDM_MZp1000_MA0300_re.root"){
    
    //inputFile = "../2HDMfullSimFile/bb2HDM_TMVA10K_MZp600_MA0300.root";
    inputFile = "../../QCDtestBGrootfile/NCUGlobalTuples_76.root";
    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    //gStyle->SetTitleAlign(33);
    //gStyle->SetTitleX(.99);
    gStyle->SetOptStat(0001111101.);
    //gStyle->SetOptStat(0);

    //get TTree from file ...
    TreeReader data(inputFile.data());
    bool isBG = false;
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2500||Zpmass.Atof()<500 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "1000";
        A0mass = "300";
        isBG = true;
    }
    vector<vector<int>> genPar;
# if doGenMatch != 0 || corr != 0
    if (!isBG) genPar = genMatch_base(inputFile.data());
#endif    
    TCanvas* c1       = new TCanvas("c1","",600,600);
    const string suf[4] = {"L","M","T","NULL"}; 
    const int njets = 10;
    const int ncom = 50;
    //const float x2MaxCut = 50;
    const float x2MaxCut = 100;
    const int   maxSaveCom = 50;
    const float CISVV2CUT_L = 0.5426;
    const float CISVV2CUT_M = 0.8484;
    const float CISVV2CUT_T = 0.9535;
    
    TH1F* h_allEvent  = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    TH1F* h_NgoodJets = new TH1F("h_NgoodJets", "h_NgoodJets", njets,-0.5,njets-0.5);
    TH1F* h_NbJets = new TH1F("h_NbJets", "h_NbJets", njets,-0.5,njets-0.5);
    TH1F* h_Ncom      = new TH1F("h_Ncom", "h_Ncom", ncom,-0.5,ncom-0.5);
    // Mass
    TH1F* h_hM        = new TH1F("h_higgsM", "h_higgsM", 70,60,200);
    TH1F* h_a0M       = new TH1F("h_a0M", "h_A0M", 50,0,500);
    TH1F* h_zpM       = new TH1F("h_ZpM", "h_ZpM", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    TH1F* h_hPt       = new TH1F("h_higgsPt", "h_higgsPt", 60,0,1200);
    TH1F* h_a0Pt      = new TH1F("h_a0Pt", "h_A0Pt", 60,0,1200);
    // study the correlation of CISVV2 and different variables
    TH2F* h_2d_hpt    = new TH2F("h_2d_hpt","h_2d_hpt",60,0,1200,100,0,1);
    TH2F* h_2d_hx2    = new TH2F("h_2d_hx2","h_2d_hx2",50,0,50,100,0,1);
    TH2F* h_2d_hptas  = new TH2F("h_2d_hptas","h_2d_hptas",40,0,2,100,0,1);
    TH2F* h_2d_hptsd  = new TH2F("h_2d_hptsd","h_2d_hptsd",35,0,0.7,100,0,1);
    TH2F* h_2d_hdeltaR  = new TH2F("h_2d_hdeltaR","h_2d_hdeltaR",40,0,4,100,0,1);
    TH2F* h_2d_hdeltaEta = new TH2F("h_2d_hdeltaEta","h_2d_hdeltaEta",40,0,4,100,0,1);
    TH2F* h_2d_hdeltaPhi = new TH2F("h_2d_hdeltaPhi","h_2d_hdeltaPhi",32,0,3.2,100,0,1);
    TH1F* h_1d_hpt[3];
    TH1F* h_1d_hx2[3];
    TH1F* h_1d_hptas[3],*h_1d_hptsd[3];
    TH1F* h_1d_hDeltaR[3],*h_1d_hDeltaEta[3], *h_1d_hDeltaPhi[3];
    string suffixA[3] = {"_fail","_pass","_ratio"};
    for (int i=0;i<3;i++) {
        string sPt = "h_1d_hpt" + suffixA[i];
        string sx2 = "h_1d_hx2" + suffixA[i];
        string sptas = "h_1d_hptas" + suffixA[i];
        string sptsd = "h_1d_hptsd" + suffixA[i];
        string sdeltaR = "h_1d_hDeltaR" + suffixA[i];
        string sdeltaEta = "h_1d_hDeltaEta" + suffixA[i];
        string sdeltaPhi = "h_1d_hDeltaPhi" + suffixA[i];

        h_1d_hpt[i]   = new TH1F(sPt.data(),sPt.data(),60,0,1200);
        h_1d_hx2[i]   = new TH1F(sx2.data(),sx2.data(),50,0,50);
        h_1d_hptas[i] = new TH1F(sptas.data(),sptas.data(),40,0,2);
        h_1d_hptsd[i] = new TH1F(sptsd.data(),sptsd.data(),35,0,0.7);
        h_1d_hDeltaR[i] = new TH1F(sdeltaR.data(),sdeltaR.data(),40,0,4);
        h_1d_hDeltaEta[i] = new TH1F(sdeltaEta.data(),sdeltaEta.data(),40,0,4);
        h_1d_hDeltaPhi[i] = new TH1F(sdeltaPhi.data(),sdeltaPhi.data(),32,0,3.2);
    }
    // Assymetry
    TH1F* h_hPtAs     = new TH1F("h_hPtAs","h_higgsPtAssymetry",40,0,2);
    TH1F* h_a0PtAs    = new TH1F("h_a0PtAs","h_A0PtAssymetry",40,0,2);
    TH1F* h_zpPtAs    = new TH1F("h_zpPtAs","h_ZpPtAssymetry",60,0,3.0);

    TH1F* h_hPtSD     = new TH1F("h_hPtSD","h_higgsPtSD",40,0,0.7);
    TH1F* h_a0PtSD    = new TH1F("h_a0PtSD","h_A0PtSD",30,0,0.6);
    TH1F* h_zpPtSD    = new TH1F("h_zpPtSD","h_ZpPtSD",30,0.1,0.7);
    // CISVV2 & weight
    TH1F *h_CISVV2_H[2], *h_CISVV2_A0[2];
    for (int i=0;i<2;i++) {
        h_CISVV2_H[i]     = new TH1F(Form("h_CISVV2_H_%d",i),Form("h_CISVV2_H_%d",i),50,0.5,1);
        h_CISVV2_A0[i]    = new TH1F(Form("h_CISVV2_A0_%d",i),Form("h_CISVV2_A0_%d",i),50,0.5,1);
        
    }
    // Delta R, theta, phi
    TH1F* h_hDeltaR         = new TH1F("h_HiggstobbDeltaR", "h_HiggstobbDeltaR", 40,0,4);
    TH1F* h_a0DeltaR        = new TH1F("h_A0tobbDeltaR", "h_A0tobbDeltaR", 40,0,4);
    TH1F* h_zpDeltaR        = new TH1F("h_ZptoHA0DeltaR", "h_ZptoHA0DeltaR_reco", 35,1,4.5);

    TH1F* h_hDeltaEta       = new TH1F("h_HiggstobbDeltaEta", "h_HiggstobbDeltaEta", 40,0,4);
    TH1F* h_a0DeltaEta      = new TH1F("h_A0tobbDeltaEta", "h_A0tobbDeltaEta", 40,0,4);
    TH1F* h_zpDeltaEta      = new TH1F("h_ZptoHA0DeltaEta", "h_ZptoHA0DeltaEta", 40,0,4);

    TH1F* h_hDeltaPhi       = new TH1F("h_HiggstobbDeltaPhi", "h_HiggstobbDeltaPhi", 32,0,3.2);
    TH1F* h_a0DeltaPhi      = new TH1F("h_A0tobbDeltaPhi", "h_A0obbDeltaPhi", 32,0,3.2);
    TH1F* h_zpDeltaPhi      = new TH1F("h_ZptoHA0DeltaPhi", "h_ZptoHA0DeltaPhi", 32,0,3.2);
    
    // declare test TH1 and TH2
    
    // test declare
    TH1F* h_Hindex1 = new TH1F("h_Hindex1","h_Hindex1",4,0.5,4.5);
    TH1F* h_Hindex2 = new TH1F("h_Hindex2","h_Hindex2",4,0.5,4.5);
    TH1F* h_A0index1 = new TH1F("h_A0index1","h_A0index1",4,0.5,4.5);
    TH1F* h_A0index2 = new TH1F("h_A0index2","h_A0index2",4,0.5,4.5);
    TH1F* h_HindexDis = new TH1F("h_HindexDis","h_HindexDis",16,-0.5,15.5);
    TH1F* h_HindexDis2 = new TH1F("h_HindexDis2","h_HindexDis",51,-0.5,50.5);
    TH1F* h_A0indexDis = new TH1F("h_A0indexDis","h_A0indexDis",16,-0.5,15.5);
    TH1F* h_A0indexDis2 = new TH1F("h_A0indexDis2","h_A0indexDis",51,-0.5,50.5);
    TH2F* h_2d_Hindex    = new TH2F("h_2d_Hindex","h_2d_Hindex",10,-0.5,9.5,15,-0.5,14.5);
    TH2F* h_2d_A0index    = new TH2F("h_2d_A0index","h_2d_A0index",10,-0.5,9.5,15,-0.5,14.5);
    TH2F* h_2d_HLMT    = new TH2F("h_2d_HLMMT","h_2d_HgenMatchLMT",4,0,4,4,0,4);
    TH2F* h_2d_A0LMT    = new TH2F("h_2d_A0LMMT","h_2d_A0genMatchLMT",4,0,4,4,0,4);
    
    TH1F* h_HmatchIndex1 = new TH1F("h_HmatchIndex1","",9,-0.5,8.5);
    TH1F* h_HmatchIndex2 = new TH1F("h_HmatchIndex2","",9,-0.5,8.5);
    TH1F* h_A0matchIndex1 = new TH1F("h_A0matchIndex1","",9,-0.5,8.5);
    TH1F* h_A0matchIndex2 = new TH1F("h_A0matchIndex2","",9,-0.5,8.5);
    TH1F* h_HmatchpassIndex1 = new TH1F("h_HmatchpassIndex1","",9,-0.5,8.5);
    TH1F* h_HmatchpassIndex2 = new TH1F("h_HmatchpassIndex2","",9,-0.5,8.5);
    TH1F* h_A0matchpassIndex1 = new TH1F("h_A0matchpassIndex1","",9,-0.5,8.5);
    TH1F* h_A0matchpassIndex2 = new TH1F("h_A0matchpassIndex3","",9,-0.5,8.5);
    TH1F* h_hbbMatchDeltaR = new TH1F("h_hbbMatchDeltaR","h_hbbMatchDeltaR",30,0,3);
    TH1F* h_hbbMatchDeltaEta = new TH1F("h_hbbMatchDeltaEta","h_hbbMatchDeltaEta",30,0,3);
    TH1F* h_hbbMatchPtAs = new TH1F("h_hbbMatchPtAs","h_hbbMatchPtAs",30,0,3);
    TH1F* h_A0bbMatchDeltaR = new TH1F("h_A0bbMatchDeltaR","h_A0bbMatchDeltaR",30,0,3);
    TH1F* h_A0bbMatchDeltaEta = new TH1F("h_A0bbMatchDeltaEta","h_A0bbMatchDeltaEta",30,0,3);
    TH1F* h_A0bbMatchPtAs = new TH1F("h_A0bbMatchPtAs","h_A0bbMatchPtAs",30,0,3);
    TH1F* h_A0HMatchDeltaR = new TH1F("h_A0HMatchDeltaR","h_A0HMatchDeltaR",30,2,5);
    TH1F* h_A0HMatchDeltaEta = new TH1F("h_A0HMatchDeltaEta","h_A0HMatchDeltaEta",30,0,3);
    TH1F* h_A0HMatchPtAs = new TH1F("h_A0HMatchPtAs","h_A0HMatchPtAs",30,0,3);
    TH1F* h_LeadPtAs = new TH1F("h_LeadPtAs","h_LeadPtAs",30,0,3);
    TH1F* h_LeadHM = new TH1F("h_LeadHM","h_LeadHM",40,80,160);
    TH1F* h_LeadA0M = new TH1F("h_LeadA0M","h_LeadA0M",50,200,400);
    // end of test
    Int_t nPass[20]={0};
    const int nSel = 1;
    Int_t nCorrectness[20]={0};
    Int_t n5Correctness[20]={0};
    Bool_t isTag, isAntiTag;
    Int_t nGoodJets;
    Float_t CISVV2_Hb1, CISVV2_Hb2, CISVV2_A0b1, CISVV2_A0b2;
    Float_t Mh, MA0;
    Float_t hPt, hDeltaR, hDeltaPhi, hDeltaEta, hptAs, hsdAs;
    Float_t a0Pt, a0DeltaR, a0DeltaPhi, a0DeltaEta, a0ptAs, a0sdAs;
    Int_t nCisT, nCisM, nCisL;
    TLorentzVector *vh1,*vh2,*va01,*va02;
    vh1 = new TLorentzVector();
    vh2 = new TLorentzVector();
    va01 = new TLorentzVector();
    va02 = new TLorentzVector();
    TFile *ftree;
    TTree *tree;
    if (saveTree) {
        ftree = new TFile(Form("treeBG_%d.root",w),"recreate");
        tree  = new TTree("tree","tree");
        tree->Branch("isTag",&isTag,"isTag/O");
        tree->Branch("isAntiTag",&isAntiTag,"isAntiTag/O");
        tree->Branch("nPassL",&nCisL,"nPassL/I");
        tree->Branch("nPassM",&nCisM,"nPassM/I");
        tree->Branch("nPassT",&nCisT,"nPassT/I");
        tree->Branch("CISVV2_Hb1",&CISVV2_Hb1,"CISVV2_Hb1/F");
        tree->Branch("CISVV2_Hb2",&CISVV2_Hb2,"CISVV2_Hb2/F");
        tree->Branch("CISVV2_A0b1",&CISVV2_A0b1,"CISVV2_A0b1/F");
        tree->Branch("CISVV2_A0b2",&CISVV2_A0b2,"CISVV2_A0b2/F");
        tree->Branch("Mh",&Mh,"Mh/F");
        tree->Branch("MA0",&MA0,"MA0/F");
        tree->Branch("hPt",&hPt,"hPt/F");
        tree->Branch("A0Pt",&a0Pt,"A0Pt/F");
        tree->Branch("hDeltaR",&hDeltaR,"hDeltaR/F");
        tree->Branch("hDeltaPhi",&hDeltaPhi,"hDeltaPhi/F");
        tree->Branch("hDeltaEta",&hDeltaEta,"hDeltaEta/F");
        tree->Branch("A0DeltaR",&a0DeltaR,"A0DeltaR/F");
        tree->Branch("A0DeltaPhi",&a0DeltaPhi,"A0DeltaPhi/F");
        tree->Branch("A0DeltaEta",&a0DeltaEta,"A0DeltaEta/F");
        tree->Branch("hptas",&hptAs,"hptAs/F");
        tree->Branch("hsdas",&hsdAs,"hsdAs/F");
        tree->Branch("A0ptas",&a0ptAs,"A0ptAs/F");
        tree->Branch("A0sdas",&a0sdAs,"A0sdAs/F");
        tree->Branch("higgs jet1","TLorentzVector",&vh1,16000);
        tree->Branch("higgs jet2","TLorentzVector",&vh2,16000);
        tree->Branch("A0 jet1","TLorentzVector",&va01,16000);
        tree->Branch("A0 jet2","TLorentzVector",&va02,16000);
    }
    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            h_2d_HLMT->Fill(suf[i].data(),suf[j].data(),0);
            h_2d_A0LMT->Fill(suf[i].data(),suf[j].data(),0);
        }
    }
    Int_t nPTTTT = 0, nPTTLL = 0, nPLLLL = 0;
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %500 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        if (jEntry>5000) break;
        data.GetEntry(jEntry);
        nPass[0]++;
        h_allEvent->Fill(1);
        isTag = false;
        isAntiTag = false;
        //0. has a good vertex
        int nGenJet =  data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[1]++;
        
        vector<bool>& vPassID_L = *(vector<bool>*) data.GetPtr("THINjetPassIDLoose");
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float *CISVV2 = data.GetPtrFloat("THINjetCISVV2");
        
        vector<TLorentzVector*> genHA0Par;
# if doGenMatch != 0 || corr !=0 
        if (!isBG) {
            genHA0Par.clear();
            TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
            for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        }
# endif
        // find good jet
        TLorentzVector* thisJet, *thatJet;
        vector<TLorentzVector*> LeadHA0Jet;
        vector<float> LeadHA0window_CISVV2;
        vector<float> LeadHA0_CISVV2;
        int matchInd[4] = {-1,-1,-1,-1};
        int matchOrderInd[4] = {-1,-1,-1,-1};
        int orderInd[4];
        TLorentzVector *HA0v4[4];
        // select 4 or 5 jet
        int nbjet = 0;
        int nCisvv2_L = 0, nCisvv2_M = 0, nCisvv2_T = 0;
        int nLjet = 0, nMjet = 0, nTjet = 0;
        int Hind[10][4];
        // loop thin Jet
        for (int i=0;i<nGenJet;i++) {
            thisJet = (TLorentzVector*)genjetP4->At(i);
            if (!isBG) {
                if (genHA0Par[0]->DeltaR(*thisJet)<0.4&&matchInd[0]<0) {matchInd[0]=i;HA0v4[0]=thisJet;}
                if (genHA0Par[1]->DeltaR(*thisJet)<0.4&&matchInd[1]<0) {matchInd[1]=i;HA0v4[1]=thisJet;}
                if (genHA0Par[2]->DeltaR(*thisJet)<0.4&&matchInd[2]<0) {matchInd[2]=i;HA0v4[2]=thisJet;}
                if (genHA0Par[3]->DeltaR(*thisJet)<0.4&&matchInd[3]<0) {matchInd[3]=i;HA0v4[3]=thisJet;}
            }
            if (thisJet->Pt()<30) continue;
            if (abs(thisJet->Eta())>2.4) continue;
            if (!vPassID_L[i]) continue;
            if (CISVV2[i]<CISVV2CUT_L) continue;
            // cal number of b jet with L, M, T
            if (CISVV2[i]>CISVV2CUT_T) {
                nCisvv2_T++;
            }
            else if (CISVV2[i]>CISVV2CUT_M) {
                nCisvv2_M++;
            }
            else if (CISVV2[i]>CISVV2CUT_L) {
                nCisvv2_L++;
            }
            // select TMML
            if (nMjet<4&&CISVV2[i]>CISVV2CUT_M) {
                Hind[0][nMjet] = i;
                nMjet++;
                LeadHA0Jet.push_back(thisJet);
                LeadHA0_CISVV2.push_back(CISVV2[i]);
            }
            nbjet++;
        
        }
        if (nMjet<2) continue;
        nPass[2]++;
        h_NbJets->Fill(nbjet);
        nCisT = nCisvv2_T;
        nCisM = nCisvv2_T + nCisvv2_M;
        nCisL = nCisvv2_T + nCisvv2_M + nCisvv2_L;
        // genmatching part
        if (matchInd[0]>=0) {h_HmatchIndex1->Fill(matchInd[0]);h_CISVV2_H[0]->Fill(CISVV2[matchInd[0]]);}
        if (matchInd[1]>=0) {h_HmatchIndex2->Fill(matchInd[1]);h_CISVV2_H[1]->Fill(CISVV2[matchInd[1]]);}
        if (matchInd[2]>=0) {h_A0matchIndex1->Fill(matchInd[2]);h_CISVV2_A0[0]->Fill(CISVV2[matchInd[2]]);}
        if (matchInd[3]>=0) {h_A0matchIndex2->Fill(matchInd[3]);h_CISVV2_A0[1]->Fill(CISVV2[matchInd[3]]);}
        if (matchInd[0]>=0&&matchInd[1]>=0) {
            h_hbbMatchDeltaR->Fill(HA0v4[0]->DeltaR(*HA0v4[1]));
            h_hbbMatchDeltaEta->Fill(abs(HA0v4[0]->Eta()-HA0v4[1]->Eta()));
            h_hbbMatchPtAs->Fill(ptAssymetry(HA0v4[0],HA0v4[1]));
            if (matchInd[0]<matchInd[1])    h_2d_Hindex->Fill(matchInd[0],matchInd[1]);
            else                            h_2d_Hindex->Fill(matchInd[1],matchInd[0]);
        }
        if (matchInd[2]>=0&&matchInd[3]>=0) {
            h_A0bbMatchDeltaR->Fill(HA0v4[2]->DeltaR(*HA0v4[3]));
            h_A0bbMatchDeltaEta->Fill(abs(HA0v4[2]->Eta()-HA0v4[3]->Eta()));
            h_A0bbMatchPtAs->Fill(ptAssymetry(HA0v4[2],HA0v4[3]));
            if (matchInd[2]<matchInd[3])    h_2d_A0index->Fill(matchInd[2],matchInd[3]);
            else                            h_2d_A0index->Fill(matchInd[3],matchInd[2]);
        
        }
        // veto event not match to genPar
        bool isMatch = matchInd[0]>=0&&matchInd[1]>=0&&matchInd[2]>=0&&matchInd[3]>=0;
        if (isMatch) {
            TLorentzVector vHA0[2];
            for (int i=0;i<2;i++) vHA0[i] = *HA0v4[2*i]+*HA0v4[2*i+1];
            h_A0HMatchDeltaR->Fill(vHA0[0].DeltaR(vHA0[1]));
            h_A0HMatchDeltaEta->Fill(abs(vHA0[0].Eta()-vHA0[1].Eta()));
            h_A0HMatchPtAs->Fill(ptAssymetry(&vHA0[0],&vHA0[1]));
            h_HindexDis->Fill(matchInd[0]*matchInd[0]+matchInd[1]*matchInd[1]); 
            h_A0indexDis->Fill(matchInd[2]*matchInd[2]+matchInd[3]*matchInd[3]); 
            h_HindexDis2->Fill(matchInd[0]*matchInd[0]+matchInd[1]*matchInd[1]); 
            h_A0indexDis2->Fill(matchInd[2]*matchInd[2]+matchInd[3]*matchInd[3]); 
            int nMatch = 0;
            for (int i=0;i<4;i++) {
                thisJet = (TLorentzVector*)genjetP4->At(i);
                if (thisJet->Pt()<30) continue;
                if (abs(thisJet->Eta())>2.4) continue;
                if (!vPassID_L[i]) continue;
                if (CISVV2[i]<CISVV2CUT_L) continue;
                nMatch++;
            }
            if (nMatch==4) {
                h_HmatchpassIndex1->Fill(matchInd[0]);
                h_HmatchpassIndex2->Fill(matchInd[1]);
                h_A0matchpassIndex1->Fill(matchInd[2]);
                h_A0matchpassIndex2->Fill(matchInd[3]);
            }
        }
        // end of genmatch
        if (LeadHA0Jet.size()==4) isTag = true;
        if (nMjet<4) isAntiTag = true;
        if (!isAntiTag&&!isTag) continue;
        nPass[3]++;
        // TMML
        if (isMatch) {
            int HId[2] = {3,3};
            int A0Id[2] = {3,3};
            for (int i=0;i<2;i++) {
                if (CISVV2[matchInd[i]]>CISVV2CUT_T) HId[i] = 2; 
                else if (CISVV2[matchInd[i]]>CISVV2CUT_M) HId[i] = 1; 
                else if (CISVV2[matchInd[i]]>CISVV2CUT_L) HId[i] = 0; 
                if (CISVV2[matchInd[i+2]]>CISVV2CUT_T) A0Id[i] = 2; 
                else if (CISVV2[matchInd[i+2]]>CISVV2CUT_M) A0Id[i] = 1; 
                else if (CISVV2[matchInd[i+2]]>CISVV2CUT_L) A0Id[i] = 0; 
            }
            if (HA0v4[0]->Pt()<HA0v4[1]->Pt()) swap(HId[0],HId[1]);
            if (HA0v4[2]->Pt()<HA0v4[3]->Pt()) swap(A0Id[0],A0Id[1]);
            h_2d_HLMT->Fill(suf[HId[0]].data(),suf[HId[1]].data(),1);
            h_2d_A0LMT->Fill(suf[A0Id[0]].data(),suf[A0Id[1]].data(),1);
        }
        
        
        float mindeMH = 999, mindeMH_range = 999, mindeMA0_5 = 999;
        // 4 b jet, H jet, A0 jet
        TLorentzVector *LeadDeltaRHA0Jet[6];
        TLorentzVector *LeadSelHA0Jet[10][6];
        for (int i=0;i<3;i++) for (int j=0;j<4;j++) Hind[i][j] = 0;
        bool inWindow = false;
        // loop for higgs and A0
        if (isAntiTag) {
            TLorentzVector hj = *LeadHA0Jet[0] + *LeadHA0Jet[1];
            if (hj.M()>135||hj.M()<115) continue;
            while(LeadHA0Jet.size()<2) {
                LeadHA0Jet.pop_back();
                LeadHA0_CISVV2.pop_back();
            }
            // search fail jet
            vector <TLorentzVector*> failJet;
            vector <int> failJetInd;
            for (int i=0;i<nGenJet;i++) {
                thisJet = (TLorentzVector*)genjetP4->At(i);
                if (thisJet->Pt()<30) continue;
                if (abs(thisJet->Eta())>2.4) continue;
                if (!vPassID_L[i]) continue;
                if (CISVV2[i]<CISVV2CUT_L) {
                    failJet.push_back(thisJet);
                    failJetInd.push_back(i);
                }
            }
            //for (int n=0;n<5;n++) for (int i=0;i<n;i++) {
            //  int j=n-i;if (i==j) continue;cout << i << j << "  "<< n << endl;}
            for (int n=0;n<failJet.size();n++) {
                for (int i=0;i<=n/2;i++) {
                    // i Leading, j trailing
                    int j = n - i;
                    if (i==Hind[0][0]||i==Hind[0][1]) continue;
                    if (j==Hind[0][1]||j==Hind[0][0]) continue;
                    if (i>=j) continue;
                    TLorentzVector vec = *failJet[i] + *failJet[j];
                    if (vec.M()>200) {
                        LeadHA0Jet.push_back(failJet[i]);
                        LeadHA0Jet.push_back(failJet[j]);
                        LeadHA0_CISVV2.push_back(CISVV2[i]);
                        LeadHA0_CISVV2.push_back(CISVV2[j]);
                        inWindow = true;
                        Hind[0][2] = i;
                        Hind[0][3] = j;
                        break;
                    }
                }
                if (inWindow) break;
            }
            if (LeadHA0Jet.size()!=4||!inWindow) continue;
        }
        nPass[4]++;
        if (isTag) {
            for (int i=0,k=0,g=0;i<4;i++) {
                for (int j=0;j<i;j++) {
                    arrangeIndex(i,j,k,g);
                    TLorentzVector hj = *LeadHA0Jet[i] + *LeadHA0Jet[j];
                    TLorentzVector a0j = *LeadHA0Jet[k] + *LeadHA0Jet[g];
                    if (hj.M()>a0j.M()) swap(hj,a0j);
                    if (a0j.M()<200) continue;
                    if (abs(hj.M()-125)<mindeMH_range&&hj.M()<135&&hj.M()>115) {
                        Hind[0][0] = i; Hind[0][1] = j;
                        Hind[0][2] = k; Hind[0][3] = g;
                        mindeMH_range = abs(hj.M()-125);
                        inWindow = true;
                    }
                }
            }
            if (!inWindow) continue;
        }
        nPass[5]++;
        cout <<jEntry << "  "<< isTag << isAntiTag << endl;
        //if (inWindow) arrangeIndex(Hind[0][0],Hind[0][1],Hind[0][2],Hind[0][3]);
        // convert index to jet
        cout << inWindow << "\t" << Hind[0][0] << " " << Hind[0][1] << " " << Hind[0][2] << " " << Hind[0][3] << endl;
        // swap h and a0
        for (int j=0;j<4;j++) LeadSelHA0Jet[0][j] = LeadHA0Jet[Hind[0][j]]; 
        LeadSelHA0Jet[0][4] = new TLorentzVector(*LeadSelHA0Jet[0][0]+*LeadSelHA0Jet[0][1]);
        LeadSelHA0Jet[0][5] = new TLorentzVector(*LeadSelHA0Jet[0][2]+*LeadSelHA0Jet[0][3]);
        if (LeadSelHA0Jet[0][4]->M()>LeadSelHA0Jet[0][5]->M()) {
            swap(LeadSelHA0Jet[0][0],LeadSelHA0Jet[0][2]);
            swap(LeadSelHA0Jet[0][1],LeadSelHA0Jet[0][3]);
            swap(LeadSelHA0Jet[0][4],LeadSelHA0Jet[0][5]);
            swap(Hind[0][0],Hind[0][2]);
            swap(Hind[0][1],Hind[0][3]);
        }
        // fill hist
        h_hM->Fill(LeadSelHA0Jet[0][4]->M());
        h_hPt->Fill(LeadSelHA0Jet[0][4]->Pt());
        h_hDeltaR->Fill(LeadSelHA0Jet[0][0]->DeltaR(*LeadSelHA0Jet[0][1]));
        h_hDeltaEta->Fill(abs(LeadSelHA0Jet[0][0]->Eta()-LeadSelHA0Jet[0][1]->Eta()));
        h_hDeltaPhi->Fill(LeadSelHA0Jet[0][0]->DeltaPhi(*LeadSelHA0Jet[0][1]));
        h_hPtAs->Fill(ptAssymetry(LeadSelHA0Jet[0][0],LeadSelHA0Jet[0][1]));
        h_hPtSD->Fill(softDropAs(LeadSelHA0Jet[0][0],LeadSelHA0Jet[0][1]));
        // end of filling
        CISVV2_Hb1 = LeadHA0_CISVV2[0];
        CISVV2_Hb2 = LeadHA0_CISVV2[1];
        CISVV2_A0b1 = LeadHA0_CISVV2[2];
        CISVV2_A0b2 = LeadHA0_CISVV2[3];
        Mh = (*LeadSelHA0Jet[0][0]+*LeadSelHA0Jet[0][1]).M();
        MA0 = (*LeadSelHA0Jet[0][2]+*LeadSelHA0Jet[0][3]).M();
        hPt = (*LeadSelHA0Jet[0][0]+*LeadSelHA0Jet[0][1]).Pt();
        a0Pt = (*LeadSelHA0Jet[0][2]+*LeadSelHA0Jet[0][3]).Pt();
        hDeltaR = LeadSelHA0Jet[0][0]->DeltaR(*LeadSelHA0Jet[0][1]);
        hDeltaPhi = LeadSelHA0Jet[0][0]->DeltaPhi(*LeadSelHA0Jet[0][1]);
        hDeltaEta = abs(LeadSelHA0Jet[0][0]->Eta()-LeadSelHA0Jet[0][1]->Eta());
        a0DeltaR = LeadSelHA0Jet[0][2]->DeltaR(*LeadSelHA0Jet[0][3]);
        a0DeltaPhi = LeadSelHA0Jet[0][2]->DeltaPhi(*LeadSelHA0Jet[0][3]);
        a0DeltaEta = abs(LeadSelHA0Jet[0][2]->Eta()-LeadSelHA0Jet[0][3]->Eta());
        hptAs = ptAssymetry(LeadSelHA0Jet[0][0],LeadSelHA0Jet[0][1]);
        hsdAs = softDropAs(LeadSelHA0Jet[0][0],LeadSelHA0Jet[0][1]);
        a0ptAs = ptAssymetry(LeadSelHA0Jet[0][2],LeadSelHA0Jet[0][3]);
        a0sdAs = softDropAs(LeadSelHA0Jet[0][2],LeadSelHA0Jet[0][3]);
        *vh1 = *LeadSelHA0Jet[0][0];
        *vh2 = *LeadSelHA0Jet[0][1];
        *va01 = *LeadSelHA0Jet[0][2];
        *va02 = *LeadSelHA0Jet[0][3];
        if (saveTree) tree->Fill();
        if (!isMatch) continue;
        nPass[6]++;
        if (inWindow) nPass[7]++;
        nCorrectness[0]++;
        if (!isBG) if (isCorrect(LeadSelHA0Jet[1],&genHA0Par[0])) nCorrectness[1]++;
        if (true||LeadSelHA0Jet[0][4]->M()>LeadSelHA0Jet[0][5]->M()) {
            delete LeadSelHA0Jet[0][4];
            delete LeadSelHA0Jet[0][5];
        }
    } // end of loop over entries
    h_2d_Hindex->GetXaxis()->SetTitle("leading H b jet"); 
    h_2d_Hindex->GetYaxis()->SetTitle("subleading H b jet"); 
    h_2d_A0index->GetXaxis()->SetTitle("leading A0 b jet"); 
    h_2d_A0index->GetYaxis()->SetTitle("subleading A0 b jet");
    h_2d_HLMT->GetXaxis()->SetTitle("leading H b jet");
    h_2d_HLMT->GetYaxis()->SetTitle("subleading H b jet");
    h_2d_A0LMT->GetXaxis()->SetTitle("leading A0 b jet");
    h_2d_A0LMT->GetYaxis()->SetTitle("subleading A0 b jet");
    h_2d_HLMT->Scale(1/h_2d_HLMT->Integral());
    h_2d_A0LMT->Scale(1/h_2d_A0LMT->Integral());
    
    float nTotal = data.GetEntriesFast();
    if (!isBG&&false) {
        h_NbJets->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf(");
        h_2d_Hindex->Draw("colz");
        c1->Print("match_MA0600_LMT.pdf");
        h_2d_A0index->Draw("colz");
        c1->Print("match_MA0600_LMT.pdf");
        h_CISVV2_H[0]->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_CISVV2_H[1]->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_CISVV2_A0[0]->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_CISVV2_A0[1]->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_HmatchIndex1->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_HmatchIndex2->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0matchIndex1->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0matchIndex2->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.16);
        h_2d_HLMT->Draw("colz text");
        c1->Print("match_MA0600_LMT.pdf");
        h_2d_A0LMT->Draw("colz text");
        c1->Print("match_MA0600_LMT.pdf");
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.1);
        h_HmatchpassIndex1->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_HmatchpassIndex2->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0matchpassIndex1->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0matchpassIndex2->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_hbbMatchDeltaR->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_hbbMatchDeltaEta->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_hbbMatchPtAs->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0bbMatchDeltaR->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0bbMatchDeltaEta->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0bbMatchPtAs->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0HMatchDeltaR->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0HMatchDeltaEta->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_A0HMatchPtAs->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_LeadPtAs->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_LeadHM->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf");
        h_LeadA0M->Draw("hist");
        c1->Print("match_MA0600_LMT.pdf)");
    }
    string nSelText[] = {"nPassSel","min DeltaR","min deMh","min PtAs","max PtAs HA0","max DeltaR HA0","higgs mass window"};
    int lastInd = -1;
    // last index is genMatching rate
    std::cout << "nTotal    = " << nTotal << std::endl;
    for(int i=0;i<20;i++) if(nPass[i]>0) {std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;lastInd = i;}
    lastInd = 4;
    efferr(nPass[lastInd],nTotal);

    string fileName; 
    if (isBG||true) {
        if (saveTree) {
            tree->Write();
            h_allEvent->Write();
            ftree->Close();
        
        }
    }
}
