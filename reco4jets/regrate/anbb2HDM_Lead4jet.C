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
#include "TMVA_regression_nu_Vali_rewrite.h"

#define saveEffi        true
#define basePtEtaCut    0
#define doGenMatch      0
#define doBTagCut       false
#define saveTree        0
#define corr            1

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
void anbb2HDM_Lead4jet(int w=0, std::string inputFile="2HDM_MZp1000_MA0300_re.root"){
    
    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    //gStyle->SetTitleAlign(33);
    //gStyle->SetTitleX(.99);
    gStyle->SetOptStat(0001111101.);
    //gStyle->SetOptStat(0);

    //get TTree from file ...
    TreeReader data(inputFile.data());
# if doGenMatch != 0 || corr != 0
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
#endif    
    bool isBG = false;
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2500||Zpmass.Atof()<500 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "1000";
        A0mass = "300";
        isBG = true;
    }
    TCanvas* c1       = new TCanvas("c1","",600,600);
    TH1F* h_allEvent  = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    
    const int njets = 10;
    const int ncom = 50;
    //const float x2MaxCut = 50;
    const float x2MaxCut = 100;
    const int   maxSaveCom = 50;
    const float CISVV2CUT_L = 0.5426;
    const float CISVV2CUT_M = 0.8484;
    const float CISVV2CUT_T = 0.9535;
    
    TH1F* h_NgoodJets = new TH1F("h_NgoodJets", "h_NgoodJets", njets,-0.5,njets-0.5);
    TH1F* h_Ncom      = new TH1F("h_Ncom", "h_Ncom", ncom,-0.5,ncom-0.5);
    // Mass
    TH1F* h_hM        = new TH1F("h_higgsM", "h_higgsM", 70,60,200);
    TH1F* h_hM_sw     = new TH1F("h_higgsM_sw", "h_higgsM_sw", 70,60,200);
    TH1F* h_hM_ws     = new TH1F("h_higgsM_ws", "h_higgsM_ws", 70,60,200);
    TH1F* h_hM_LT     = new TH1F("h_higgsM_ws_LT", "h_higgsM_ws_LT", 70,60,200);
    TH1F* h_a0M       = new TH1F("h_a0M", "h_A0M", 50,0,500);
    TH1F* h_a0M_sw    = new TH1F("h_A0M_sw", "h_A0M_sw", 50,0,500);
    TH1F* h_a0M_ws    = new TH1F("h_A0M_ws", "h_A0M_ws", 50,0,500);
    TH1F* h_a0M_LT    = new TH1F("h_A0M_ws_LT", "h_A0M_ws_LT", 50,0,500);
    TH1F* h_zpM       = new TH1F("h_ZpM", "h_ZpM", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    TH1F* h_zpM_sw    = new TH1F("h_ZpM_sw", "h_ZpM_sw", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    TH1F* h_zpM_ws    = new TH1F("h_ZpM_ws", "h_ZpM_ws", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    TH1F* h_zpM_LT    = new TH1F("h_ZpM_ws_LT", "h_ZpM_ws_LT", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    
    TH1F* h_hPt       = new TH1F("h_higgsPt", "h_higgsPt", 60,0,1200);
    TH1F* h_a0Pt      = new TH1F("h_a0Pt", "h_A0Pt", 60,0,1200);
    TH1F* h_hPt_sw    = new TH1F("h_higgsPt_sw", "h_higgsPt_sw", 60,0,1200);
    TH1F* h_a0Pt_sw   = new TH1F("h_a0Pt_sw", "h_A0Pt_sw", 60,0,1200);
    TH1F* h_hPt_ws    = new TH1F("h_higgsPt_ws", "h_higgsPt_ws", 60,0,1200);
    TH1F* h_a0Pt_ws   = new TH1F("h_a0Pt_ws", "h_A0Pt_ws", 60,0,1200);
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
    TH1F* h_hPtAs_sw  = new TH1F("h_hPtAs_sw","h_higgsPtAssymetry_sw",40,0,2);
    TH1F* h_hPtAs_ws  = new TH1F("h_hPtAs_ws","h_higgsPtAssymetry_ws",40,0,2);
    TH1F* h_a0PtAs    = new TH1F("h_a0PtAs","h_A0PtAssymetry",40,0,2);
    TH1F* h_a0PtAs_sw = new TH1F("h_a0PtAs_sw","h_A0PtAssymetry_sw",40,0,2);
    TH1F* h_a0PtAs_ws = new TH1F("h_a0PtAs_ws","h_A0PtAssymetry_ws",40,0,2);
    TH1F* h_zpPtAs    = new TH1F("h_zpPtAs","h_ZpPtAssymetry",60,0,3.0);
    TH1F* h_zpPtAs_sw = new TH1F("h_zpPtAs_sw","h_ZpPtAssymetry_sw",60,0,3.0);
    TH1F* h_zpPtAs_ws = new TH1F("h_zpPtAs_ws","h_ZpPtAssymetry_ws",60,0,3.0);

    TH1F* h_hPtSD     = new TH1F("h_hPtSD","h_higgsPtSD",40,0,0.7);
    TH1F* h_hPtSD_sw  = new TH1F("h_hPtSD_sw","h_higgsPtSD_sw",40,0,0.7);
    TH1F* h_hPtSD_ws  = new TH1F("h_hPtSD_ws","h_higgsPtSD_ws",40,0,0.7);
    TH1F* h_a0PtSD    = new TH1F("h_a0PtSD","h_A0PtSD",30,0,0.6);
    TH1F* h_a0PtSD_sw = new TH1F("h_a0PtSD_sw","h_A0PtSD_sw",30,0,0.6);
    TH1F* h_a0PtSD_ws = new TH1F("h_a0PtSD_ws","h_A0PtSD_ws",30,0,0.6);
    TH1F* h_zpPtSD    = new TH1F("h_zpPtSD","h_ZpPtSD",30,0.1,0.7);
    TH1F* h_zpPtSD_sw = new TH1F("h_zpPtSD_sw","h_ZpPtSD_sw",30,0.1,0.7);
    TH1F* h_zpPtSD_ws = new TH1F("h_zpPtSD_ws","h_ZpPtSD_ws",30,0.1,0.7);
    // CISVV2 & weight
    TH1F *h_CISVV2_H[2], *h_CISVV2_A0[2];
    TH1F *h_TMVAweight_sw[2];
    TH1F *h_TMVAweight_ws[2];
    for (int i=0;i<2;i++) {
        h_CISVV2_H[i]     = new TH1F(Form("h_CISVV2_H_%d",i),Form("h_CISVV2_H_%d",i),50,0.5,1);
        h_CISVV2_A0[i]    = new TH1F(Form("h_CISVV2_A0_%d",i),Form("h_CISVV2_A0_%d",i),50,0.5,1);
        h_TMVAweight_ws[i]  = new TH1F(Form("h_TMVAweightWS_%d",i),Form("h_TMVAweightWS_%d",i),40,0.0,2);
        h_TMVAweight_sw[i]  = new TH1F(Form("h_TMVAweightSW_%d",i),Form("h_TMVAweightSW_%d",i),40,0.0,2);
        
    }
    // pf Ratio
    TH1F* h_hM_pfRation[3][3]; // high Pt
    TH1F* h_hM_pfRationMin[3][3]; 
    TH1F* h_hM_pfRationMean[3][3]; 
    TH1F* h_a0M_pfRation[3][3]; 
    TH1F* h_zpM_pfRation[3][3];
    string suffixM[3] = {"higgsM","A0M","ZpM"}; 
    string suffixB[3] = {"_L","_M","_T"};
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            string name[3], name_min[3], name_mean[3];
            for (int k=0;k<3;k++) name[k] = "h_highPt"+suffixM[k]+suffixA[i]+suffixB[j]; 
            for (int k=0;k<3;k++) name_min[k] = "h_min"+suffixM[k]+suffixA[i]+suffixB[j]; 
            for (int k=0;k<3;k++) name_mean[k] = "h_mean"+suffixM[k]+suffixA[i]+suffixB[j]; 
            h_hM_pfRationMin[i][j] = new TH1F(name_min[0].data(),name_min[0].data(),32,40,200);
            h_hM_pfRationMean[i][j] = new TH1F(name_mean[0].data(),name_mean[0].data(),32,40,200);
            h_hM_pfRation[i][j] = new TH1F(name[0].data(),name[0].data(),32,40,200);
            h_a0M_pfRation[i][j] = new TH1F(name[1].data(),name[1].data(),60,0,600);
            h_zpM_pfRation[i][j] = new TH1F(name[2].data(),name[2].data(),80,Zpmass.Atof()-400,Zpmass.Atof()+400);
            h_hM_pfRationMin[i][j]->Sumw2();
            h_hM_pfRationMean[i][j]->Sumw2();
            h_hM_pfRation[i][j]->Sumw2();
            h_a0M_pfRation[i][j]->Sumw2();
            h_zpM_pfRation[i][j]->Sumw2();
        }
    }
    // Delta R, theta, phi
    TH1F* h_hDeltaR         = new TH1F("h_HiggstobbDeltaR", "h_HiggstobbDeltaR", 40,0,4);
    TH1F* h_hDeltaR_sw      = new TH1F("h_HiggstobbDeltaR_sw", "h_HiggstobbDeltaR_sw", 40,0,4);
    TH1F* h_hDeltaR_ws      = new TH1F("h_HiggstobbDeltaR_ws", "h_HiggstobbDeltaR_ws", 40,0,4);
    TH1F* h_a0DeltaR        = new TH1F("h_A0tobbDeltaR", "h_A0tobbDeltaR", 40,0,4);
    TH1F* h_a0DeltaR_sw     = new TH1F("h_A0tobbDeltaR_sw", "h_A0tobbDeltaR_sw", 40,0,4);
    TH1F* h_a0DeltaR_ws     = new TH1F("h_A0tobbDeltaR_ws", "h_A0tobbDeltaR_ws", 40,0,4);
    TH1F* h_zpDeltaR        = new TH1F("h_ZptoHA0DeltaR", "h_ZptoHA0DeltaR_reco", 35,1,4.5);
    TH1F* h_zpDeltaR_sw     = new TH1F("h_ZptoHA0DeltaR_sw", "h_ZptoHA0DeltaR_sw", 35,1,4.5);
    TH1F* h_zpDeltaR_ws     = new TH1F("h_ZptoHA0DeltaR_ws", "h_ZptoHA0DeltaR_ws", 35,1,4.5);

    TH1F* h_hDeltaEta       = new TH1F("h_HiggstobbDeltaEta", "h_HiggstobbDeltaEta", 40,0,4);
    TH1F* h_hDeltaEta_sw    = new TH1F("h_HiggstobbDeltaEta_sw", "h_HiggstobbDeltaEta_sw", 40,0,4);
    TH1F* h_hDeltaEta_ws    = new TH1F("h_HiggstobbDeltaEta_ws", "h_HiggstobbDeltaEta_ws", 40,0,4);
    TH1F* h_a0DeltaEta      = new TH1F("h_A0tobbDeltaEta", "h_A0tobbDeltaEta", 40,0,4);
    TH1F* h_a0DeltaEta_sw   = new TH1F("h_A0tobbDeltaEta_sw", "h_A0tobbDeltaEta_sw", 40,0,4);
    TH1F* h_a0DeltaEta_ws   = new TH1F("h_A0tobbDeltaEta_ws", "h_A0tobbDeltaEta_ws", 40,0,4);
    TH1F* h_zpDeltaEta      = new TH1F("h_ZptoHA0DeltaEta", "h_ZptoHA0DeltaEta", 40,0,4);
    TH1F* h_zpDeltaEta_sw   = new TH1F("h_ZptoHA0DeltaEta_sw", "h_ZptoHA0DeltaEta_sw", 40,0,4);
    TH1F* h_zpDeltaEta_ws   = new TH1F("h_ZptoHA0DeltaEta_ws", "h_ZptoHA0DeltaEta_ws", 40,0,4);

    TH1F* h_hDeltaPhi       = new TH1F("h_HiggstobbDeltaPhi", "h_HiggstobbDeltaPhi", 32,0,3.2);
    TH1F* h_hDeltaPhi_sw    = new TH1F("h_HiggstobbDeltaPhi_sw", "h_HiggstobbDeltaPhi_sw", 32,0,3.2);
    TH1F* h_hDeltaPhi_ws    = new TH1F("h_HiggstobbDeltaPhi_ws", "h_HiggstobbDeltaPhi_ws", 32,0,3.2);
    TH1F* h_a0DeltaPhi      = new TH1F("h_A0tobbDeltaPhi", "h_A0obbDeltaPhi", 32,0,3.2);
    TH1F* h_a0DeltaPhi_sw   = new TH1F("h_A0tobbDeltaPhi_sw", "h_A0obbDeltaPhi_sw", 32,0,3.2);
    TH1F* h_a0DeltaPhi_ws   = new TH1F("h_A0tobbDeltaPhi_ws", "h_A0obbDeltaPhi_ws", 32,0,3.2);
    TH1F* h_zpDeltaPhi      = new TH1F("h_ZptoHA0DeltaPhi", "h_ZptoHA0DeltaPhi", 32,0,3.2);
    TH1F* h_zpDeltaPhi_sw   = new TH1F("h_ZptoHA0DeltaPhi_sw", "h_ZptoHA0DeltaPhi_sw", 32,0,3.2);
    TH1F* h_zpDeltaPhi_ws   = new TH1F("h_ZptoHA0DeltaPhi_ws", "h_ZptoHA0DeltaPhi_ws", 32,0,3.2);
    
    // declare test TH1 and TH2
    
    // test declare
    TH1F* h_nComH = new TH1F("h_nComH","h_nComH",5,-0.5,4.5);
    TH1F* h_nComH_m = new TH1F("h_nComH_m","h_nComH_m",5,-0.5,4.5);
    TH1F* h_Hindex1 = new TH1F("h_Hindex1","h_Hindex1",4,0.5,4.5);
    TH1F* h_Hindex2 = new TH1F("h_Hindex2","h_Hindex2",4,0.5,4.5);
    TH1F* h_A0index1 = new TH1F("h_A0index1","h_A0index1",4,0.5,4.5);
    TH1F* h_A0index2 = new TH1F("h_A0index2","h_A0index2",4,0.5,4.5);
    TH1F* h_HindexDis = new TH1F("h_HindexDis","h_HindexDis",16,-0.5,15.5);
    TH1F* h_A0indexDis = new TH1F("h_A0indexDis","h_A0indexDis",16,-0.5,15.5);
    
    TH1F* h_HmatchIndex1 = new TH1F("h_HmatchIndex1","",9,-0.5,8.5);
    TH1F* h_HmatchIndex2 = new TH1F("h_HmatchIndex2","",9,-0.5,8.5);
    TH1F* h_A0matchIndex1 = new TH1F("h_A0matchIndex1","",9,-0.5,8.5);
    TH1F* h_A0matchIndex2 = new TH1F("h_A0matchIndex2","",9,-0.5,8.5);
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
    const int nSel = 5;
    Int_t nCorrectness[20]={0};
    Int_t n5Correctness[20]={0};
    TFile *ftree;
    TTree *tree;
    Int_t nCom, nGoodJets;
    Float_t CISVV2_Hb1[maxSaveCom], CISVV2_Hb2[maxSaveCom], CISVV2_A0b1[maxSaveCom], CISVV2_A0b2[maxSaveCom];
    Float_t TMVAweight_Hb1[maxSaveCom], TMVAweight_Hb2[maxSaveCom], TMVAweight_A0b1[maxSaveCom], TMVAweight_A0b2[maxSaveCom];
    Float_t Mh[maxSaveCom], Mh_weightLT[maxSaveCom], Mh_weight[maxSaveCom];
    Float_t hPt[maxSaveCom], hDeltaR[maxSaveCom], hDeltaPhi[maxSaveCom], hDeltaEta[maxSaveCom], hptAs[maxSaveCom], hsdAs[maxSaveCom];
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %500 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        nPass[0]++;
        h_allEvent->Fill(1);
        
        //0. has a good vertex
        int nGenJet =  data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[1]++;
        
        vector<bool>& vPassID_L = *(vector<bool>*) data.GetPtr("THINjetPassIDLoose");
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float *CISVV2 = data.GetPtrFloat("THINjetCISVV2");
        
        static vector<TLorentzVector*> genHA0Par;
# if doGenMatch != 0 || corr !=0 
        if (!isBG) {
            genHA0Par.clear();
            TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
            for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        }
# endif
        // find good jet
        static vector<int> goodJet;
        goodJet.clear();
        static TLorentzVector* thisJet, *thatJet;
        vector<TLorentzVector*> LeadHA0Jet;
        for (int i=0;i<nGenJet;i++) {
            if (!doGenMatch) {
                //if (!vPassID_L[i]) continue;
                //if (CISVV2[i]<CISVV2CUT_L) continue;
                //if (CISVV2[i]<0) continue;
                thisJet = (TLorentzVector*)genjetP4->At(i);
                //if (thisJet->Pt()<30) continue;
                //if (abs(thisJet->Eta())>2.4) continue;
            }
            goodJet.push_back(i);
            //if (goodJet.size()<4) goodJet.push_back(i);
        }
        h_NgoodJets->Fill(goodJet.size());
        int matchInd[4] = {-1,-1,-1,-1};
        int matchInd_sort[4];
        int matchOrderInd[4] = {-1,-1,-1,-1};
        int orderInd[4];
        TLorentzVector *HA0v4[4];
        // select 4 or 5 jet
        for (int i=0;i<nGenJet;i++) {
            if (genHA0Par[0]->DeltaR(*thisJet)<0.4&&matchInd[0]<0) {matchInd[0]=i;HA0v4[0]=thisJet;}
            if (genHA0Par[1]->DeltaR(*thisJet)<0.4&&matchInd[1]<0) {matchInd[1]=i;HA0v4[1]=thisJet;}
            if (genHA0Par[2]->DeltaR(*thisJet)<0.4&&matchInd[2]<0) {matchInd[2]=i;HA0v4[2]=thisJet;}
            if (genHA0Par[3]->DeltaR(*thisJet)<0.4&&matchInd[3]<0) {matchInd[3]=i;HA0v4[3]=thisJet;}
            thisJet = (TLorentzVector*)genjetP4->At(i);
            if (thisJet->Pt()<30) continue;
            if (abs(thisJet->Eta())>2.4) continue;
            if (!vPassID_L[i]) continue;
            if (CISVV2[i]<CISVV2CUT_L) continue;
            if (LeadHA0Jet.size()<5) LeadHA0Jet.push_back(thisJet);
        
        }
            
        if (matchInd[0]>=0) {h_HmatchIndex1->Fill(matchInd[0]);h_CISVV2_H[0]->Fill(CISVV2[matchInd[0]]);}
        if (matchInd[1]>=0) {h_HmatchIndex2->Fill(matchInd[1]);h_CISVV2_H[1]->Fill(CISVV2[matchInd[1]]);}
        if (matchInd[2]>=0) {h_A0matchIndex1->Fill(matchInd[2]);h_CISVV2_A0[0]->Fill(CISVV2[matchInd[2]]);}
        if (matchInd[3]>=0) {h_A0matchIndex2->Fill(matchInd[3]);h_CISVV2_A0[1]->Fill(CISVV2[matchInd[3]]);}
        if (matchInd[0]>=0&&matchInd[1]>=0) {
            h_hbbMatchDeltaR->Fill(HA0v4[0]->DeltaR(*HA0v4[1]));
            h_hbbMatchDeltaEta->Fill(abs(HA0v4[0]->Eta()-HA0v4[1]->Eta()));
            h_hbbMatchPtAs->Fill(ptAssymetry(HA0v4[0],HA0v4[1]));
        }
        if (matchInd[2]>=0&&matchInd[3]>=0) {
            h_A0bbMatchDeltaR->Fill(HA0v4[2]->DeltaR(*HA0v4[3]));
            h_A0bbMatchDeltaEta->Fill(abs(HA0v4[2]->Eta()-HA0v4[3]->Eta()));
            h_A0bbMatchPtAs->Fill(ptAssymetry(HA0v4[2],HA0v4[3]));
        
        }
        // veto event not match to genPar
        bool isMatch = false;
        if (matchInd[0]>=0&&matchInd[1]>=0&&matchInd[2]>=0&&matchInd[3]>=0) {
            isMatch = true;
            TLorentzVector vHA0[2];
            for (int i=0;i<2;i++) vHA0[i] = *HA0v4[2*i]+*HA0v4[2*i+1];
            h_A0HMatchDeltaR->Fill(vHA0[0].DeltaR(vHA0[1]));
            h_A0HMatchDeltaEta->Fill(abs(vHA0[0].Eta()-vHA0[1].Eta()));
            h_A0HMatchPtAs->Fill(ptAssymetry(&vHA0[0],&vHA0[1]));
            h_HindexDis->Fill(matchInd[0]*matchInd[0]+matchInd[1]*matchInd[1]); 
            h_A0indexDis->Fill(matchInd[2]*matchInd[2]+matchInd[3]*matchInd[3]); 

        }
        if (LeadHA0Jet.size()<4) continue;
        nPass[2]++;
        if (!isMatch) continue;
        nPass[3]++;
        // reuse matchInd
        float maxDeltaR = -1, maxDeltaR_5 = -1; 
        float maxPtAs = -1, maxPtAs_5 = -1;  
        float minPtAs = 99, minPtAs_51 = 99, minPtAs_52 = 99;
        float minDeltaR = 99, minDeltaR_51 = 99, minDeltaR_52 = 99; 
        float mindeMH = 999, mindeMH_5 = 999, mindeMA0_5 = 999;
        // 4 b jet, H jet, A0 jet
        TLorentzVector *LeadDeltaRHA0Jet[6];
        TLorentzVector *LeadSelHA0Jet[10][6];
        TLorentzVector *LeadSelHA0Jet_5[10][6];
        int Hind[10][4];
        int Hind_5[10][4];
        // loop for 5 jet H
        for (int i=0;i<5;i++) {
            if (LeadHA0Jet.size()<5 && i==4) continue;
           for (int j=0;j<i;j++) {
                if (LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j])<minDeltaR_51) {
                    minDeltaR_51 = LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j]);
                    Hind_5[0][2] = i; Hind_5[0][3] = j;
                }
                if (abs((*LeadHA0Jet[i]+*LeadHA0Jet[j]).M()-125)<mindeMH_5) {
                    mindeMH_5 = abs((*LeadHA0Jet[i]+*LeadHA0Jet[j]).M()-125);
                    Hind_5[1][0] = i; Hind_5[1][1] = j;
                }
                if (ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j])<minPtAs_51) {
                    minPtAs_51 = ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j]);    
                    Hind_5[2][2] = i; Hind_5[2][3] = j;
                }

               // 4 jet
               if (i==4) continue; 
               if (LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j])<minDeltaR) {
                    minDeltaR = LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j]);
                    Hind[0][0] = i; Hind[0][1] = j;
                }
                if (abs((*LeadHA0Jet[i]+*LeadHA0Jet[j]).M()-125)<mindeMH) {
                    mindeMH = abs((*LeadHA0Jet[i]+*LeadHA0Jet[j]).M()-125);
                    Hind[1][0] = i; Hind[1][1] = j;
                }
                if (ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j])<minPtAs) {
                    minPtAs = ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j]);    
                    Hind[2][0] = i; Hind[2][1] = j;
                }

           }
        }
        // loop for 5 jet A0
        for (int i=0;i<5;i++) {
            if (LeadHA0Jet.size()<5) break;
            for (int j=0;j<i;j++) {
                if (!isOverlap(Hind_5[0][2],Hind_5[0][3],i,j)&&LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j])<minDeltaR_52&&LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j])>minDeltaR_51) {
                    Hind_5[0][0] = i; Hind_5[0][1] = j;
                    minDeltaR_52 = LeadHA0Jet[i]->DeltaR(*LeadHA0Jet[j]);
                }
                if (!isOverlap(Hind_5[1][0],Hind_5[1][1],i,j)&&abs((*LeadHA0Jet[i]+*LeadHA0Jet[j]).M()-A0mass.Atof())<mindeMA0_5) {
                    Hind_5[1][2] = i; Hind_5[1][3] = j;
                    mindeMA0_5 = abs((*LeadHA0Jet[i]+*LeadHA0Jet[j]).M()-A0mass.Atof());
                }
                if (!isOverlap(Hind_5[2][2],Hind_5[2][3],i,j)&&ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j])<minPtAs_52&&ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j])>minPtAs_51) {
                    minPtAs_52 = ptAssymetry(LeadHA0Jet[i],LeadHA0Jet[j]);    
                    Hind_5[2][0] = i; Hind_5[2][1] = j;
                }
                TLorentzVector v1 = *LeadHA0Jet[j] + *LeadHA0Jet[i];
                for (int m=0;m<j;m++) {
                    for (int n=0;n<m;n++) {
                        if (isOverlap(i,j,n,m)) continue;
                        TLorentzVector v2 = *LeadHA0Jet[n] + *LeadHA0Jet[m];
                        if (ptAssymetry(&v1,&v2)>maxPtAs_5) {
                            maxPtAs_5 = ptAssymetry(&v1,&v2);
                            Hind_5[3][0] = i; Hind_5[3][1] = j;
                            Hind_5[3][2] = m; Hind_5[3][3] = n;
                        }
                        if (v1.DeltaR(v2)>maxDeltaR_5) {
                            maxDeltaR_5 = v1.DeltaR(v2);
                            Hind_5[4][0] = i; Hind_5[4][1] = j;
                            Hind_5[4][2] = m; Hind_5[4][3] = n;
                        }
                    }
                }
            }
        }
        // loop for 4 jet 
        for (int i=1;i<4;i++) {
           // v1 and v2 are H and A0 jet
           TLorentzVector v1 = *LeadHA0Jet[0] + *LeadHA0Jet[i];
           TLorentzVector v2;
           for (int j=1;j<4;j++) {
                if (j==i) continue;
                v2+=*LeadHA0Jet[j];
           }
           if (ptAssymetry(&v1,&v2)>maxPtAs) {
                maxPtAs = ptAssymetry(&v1,&v2);
                Hind[3][0] = 0; Hind[3][1] = i;
           }
           if (v1.DeltaR(v2)>maxDeltaR) {
                maxDeltaR = v1.DeltaR(v2);
                Hind[4][0] = 0; Hind[4][1] = i;
           }
        }
        //for (int i=0;i<nSel;i++) cout << i << "   "<< Hind_5[i][0]<<" "<<Hind_5[i][1]<<" "<<Hind_5[i][2]<<Hind_5[i][3] << endl;
        for (int i=0;i<nSel;i++) arrangeIndex(Hind[i][0],Hind[i][1],Hind[i][2],Hind[i][3]);
        // convert index to jet
        for (int i=0;i<nSel;i++) for (int j=0;j<4;j++) LeadSelHA0Jet[i][j] = LeadHA0Jet[Hind[i][j]]; 
        if (LeadHA0Jet.size()==5) for (int i=0;i<nSel;i++) for (int j=0;j<4;j++) LeadSelHA0Jet_5[i][j] = LeadHA0Jet[Hind_5[i][j]]; 
        // swap h and a0
        for (int i=0;i<nSel;i++) {
            // 4 jet
            LeadSelHA0Jet[i][4] = new TLorentzVector(*LeadSelHA0Jet[i][0]+*LeadSelHA0Jet[i][1]);
            LeadSelHA0Jet[i][5] = new TLorentzVector(*LeadSelHA0Jet[i][2]+*LeadSelHA0Jet[i][3]);
            if (LeadSelHA0Jet[i][4]->M()>LeadSelHA0Jet[i][5]->M()) {
                swap(LeadSelHA0Jet[i][0],LeadSelHA0Jet[i][2]);
                swap(LeadSelHA0Jet[i][1],LeadSelHA0Jet[i][3]);
                swap(LeadSelHA0Jet[i][4],LeadSelHA0Jet[i][5]);
                swap(Hind[i][0],Hind[i][2]);
                swap(Hind[i][1],Hind[i][3]);
            }
            // 5 jet
            if (LeadHA0Jet.size()==5) {
                LeadSelHA0Jet_5[i][4] = new TLorentzVector(*LeadSelHA0Jet_5[i][0]+*LeadSelHA0Jet_5[i][1]);
                LeadSelHA0Jet_5[i][5] = new TLorentzVector(*LeadSelHA0Jet_5[i][2]+*LeadSelHA0Jet_5[i][3]);
                if (LeadSelHA0Jet_5[i][4]->M()>LeadSelHA0Jet_5[i][5]->M()) {
                    swap(LeadSelHA0Jet_5[i][0],LeadSelHA0Jet_5[i][2]);
                    swap(LeadSelHA0Jet_5[i][1],LeadSelHA0Jet_5[i][3]);
                    swap(LeadSelHA0Jet_5[i][4],LeadSelHA0Jet_5[i][5]);
                    swap(Hind_5[i][0],Hind_5[i][2]);
                    swap(Hind_5[i][1],Hind_5[i][3]);
                    
                }
            }
        }
        nCorrectness[0]++;
        n5Correctness[0]++;
        for (int i=0;i<nSel;i++) if (isCorrect(LeadSelHA0Jet[i],&genHA0Par[0])) nCorrectness[i+1]++;
        if (LeadHA0Jet.size()==5) for (int i=0;i<nSel;i++) if (isCorrect(LeadSelHA0Jet_5[i],&genHA0Par[0])) n5Correctness[i+1]++;
        for (int i=0;i<nSel;i++) {
            delete LeadSelHA0Jet[i][4];
            delete LeadSelHA0Jet[i][5];
            if (LeadHA0Jet.size()==5) delete LeadSelHA0Jet_5[i][4];
            if (LeadHA0Jet.size()==5) delete LeadSelHA0Jet_5[i][5];
        }
    } // end of loop over entries
    
    /*
    h_zpPtSD->SetMaximum(1000);
    h_zpPtAs_ori->SetMaximum(150); 
    */
    if (saveTree) {
        tree->Write();
        h_allEvent->Write();
        ftree->Close();
    
    }
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    h_HindexDis->Draw("hist");
    c1->Print("matchIndex.pdf(");
    h_CISVV2_H[0]->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_CISVV2_H[1]->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_CISVV2_A0[0]->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_CISVV2_A0[1]->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0indexDis->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_HmatchIndex1->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_HmatchIndex2->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0matchIndex1->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0matchIndex2->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_hbbMatchDeltaR->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_hbbMatchDeltaEta->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_hbbMatchPtAs->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0bbMatchDeltaR->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0bbMatchDeltaEta->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0bbMatchPtAs->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0HMatchDeltaR->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0HMatchDeltaEta->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_A0HMatchPtAs->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_LeadPtAs->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_LeadHM->Draw("hist");
    c1->Print("matchIndex.pdf");
    h_LeadA0M->Draw("hist");
    c1->Print("matchIndex.pdf)");
    
    int lastInd = -1;
    for(int i=0;i<20;i++) if(nPass[i]>0) {std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;lastInd = i;}
    efferr(nPass[lastInd],nTotal);
    nCorrectness[nSel+1] = nPass[2]; 
    n5Correctness[nSel+1] = nPass[2];
    cout << "4 jets \t5 jets" << endl;
    for (int i=0;i<nSel+1;i++) cout << (float)nCorrectness[i]/nPass[lastInd] << "\t" << (float)n5Correctness[i]/nPass[lastInd]<<  endl;
    savenPass(nCorrectness,n5Correctness,Form("Correctness/corr_MZp%d_MA0%d.txt",Zpmass.Atoi(),A0mass.Atoi()));
}

