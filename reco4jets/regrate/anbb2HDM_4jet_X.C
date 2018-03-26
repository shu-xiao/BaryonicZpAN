#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include <TLegend.h>
#include "string"
#include "../setNCUStyle.C"
#include "../../gen2HDMsample/genMatch.C"
#include "TMVA_regression_nu_Vali_st.h"
//#include "TMVA_regression_nu_Vali_f.h"

#define saveEffi          true
#define basePtEtaCut    true
#define doGenMatch      false
#define doBTagCut       false
#define CISVV2Cut       0.5426

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
    else dePhi = abs(phi1-phi2);
    return dePhi;
}
bool sortListbyPt(vector<float> a, vector<float> b) {return a[2]>b[2];}
bool sortListbyPtZp(vector<float> a, vector<float> b) {return a[4]>b[4];}
bool sortListbyxSquare(HA0JetInfo s1, HA0JetInfo s2) {return s1.xSqure()<s2.xSqure();}
bool sortListbyWeightxSquare(HA0JetInfo s1, HA0JetInfo s2) {return s1.weightJets()->xSqure()<s2.weightJets()->xSqure();}
inline float ptAssymetry(TLorentzVector* j1, TLorentzVector* j2) {
    float minPt = (j1->Pt()>j2->Pt())? j2->Pt():j1->Pt();
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

void drawCISVV2(TH1F* h1, TH1F* h2,TH1F *h3,TH1F *h4, float weight = 1.1, int i = kBlue) {
    static float max;
    if ((h1->GetMaximum()>h2->GetMaximum())&&(h1->GetMaximum()>h3->GetMaximum())&&(h1->GetMaximum()>h4->GetMaximum())) max = h1->GetMaximum();
    else if ((h2->GetMaximum()>h1->GetMaximum())&&(h2->GetMaximum()>h3->GetMaximum())&&(h2->GetMaximum()>h4->GetMaximum())) max = h2->GetMaximum();
    else if ((h3->GetMaximum()>h2->GetMaximum())&&(h3->GetMaximum()>h1->GetMaximum())&&(h3->GetMaximum()>h4->GetMaximum())) max = h2->GetMaximum();
    else max = h4->GetMaximum();
    h1->SetMaximum(max*weight);
    h2->SetMaximum(max*weight);
    h3->SetMaximum(max*weight);
    h4->SetMaximum(max*weight);
    h3->SetLineColor(i);
    h4->SetLineColor(i);
    h1->Draw("hist");
    h2->Draw("histsame");
    h3->Draw("histsame");
    h4->Draw("histsame");
}
void draw(TH1F* h1, TH1F* h2, float weight = 1.1, int i = kBlue) {
    h1->SetMaximum((h1->GetMaximum()>h2->GetMaximum())? h1->GetMaximum()*weight : h2->GetMaximum()*weight );
    h2->SetMaximum((h1->GetMaximum()>h2->GetMaximum())? h1->GetMaximum()*weight : h2->GetMaximum()*weight );
    h1->SetLineColor(1);
    h2->SetLineColor(i);
    h1->Draw("hist");
    h2->Draw("histsame");
}
void draw(TH1F* h1, TH1F* h2, TH1F* h3, float weight = 1.1, int i = kBlue-5, int j = kBlue){
    static float max;
    if ((h1->GetMaximum()>h2->GetMaximum())&&(h1->GetMaximum()>h3->GetMaximum())) max = h1->GetMaximum();
    else if ((h2->GetMaximum()>h1->GetMaximum())&&(h2->GetMaximum()>h3->GetMaximum())) max = h2->GetMaximum();
    else max = h3->GetMaximum();
    h1->SetMaximum(max*weight);
    h2->SetMaximum(max*weight);
    h3->SetMaximum(max*weight);
    h1->SetLineColor(1);
    h2->SetLineColor(i);
    h3->SetLineColor(j);
    h1->Draw("hist");
    h2->Draw("histsame");
    h3->Draw("histsame");
}
using namespace std;
void anbb2HDM_4jet_X(int w=0, std::string inputFile="2HDM_MZp1000_MA0300_re.root"){
    
    setNCUStyle(true);
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0001111101.);
    //gStyle->SetOptStat(0);

    //get TTree from file ...
    TreeReader data(inputFile.data());
    vector<vector<int>> genPar = genMatch_base(inputFile.data());
    
    bool isBG = false;
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2500||Zpmass.Atof()<500 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "800";
        A0mass = "300";
        isBG = true;
    }
    TCanvas* c1       = new TCanvas("c1","",600,600);
    TH1F* h_allEvent  = new TH1F("h_allEvent","h_allEvent",10,-0.5,9);
    
    const int njets = 10;
    const int ncom = 50;

    TH1F* h_NgoodJets = new TH1F("h_NgoodJets", "h_NgoodJets", njets,-0.5,njets-0.5);
    TH1F* h_Ncom      = new TH1F("h_Ncom", "h_Ncom", ncom,-0.5,ncom-0.5);
    // Mass
    TH1F* h_hM        = new TH1F("h_higgsM", "h_higgsM", 70,60,200);
    TH1F* h_hM_sw     = new TH1F("h_higgsM_sw", "h_higgsM_sw", 70,60,200);
    TH1F* h_hM_ws     = new TH1F("h_higgsM_ws", "h_higgsM_ws", 70,60,200);
    TH1F* h_a0M       = new TH1F("h_a0M", "h_A0M", 50,0,500);
    TH1F* h_a0M_sw    = new TH1F("h_A0M_sw", "h_A0M_sw", 50,0,500);
    TH1F* h_a0M_ws    = new TH1F("h_A0M_ws", "h_A0M_ws", 50,0,500);
    TH1F* h_zpM       = new TH1F("h_ZpM", "h_ZpM", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    TH1F* h_zpM_sw    = new TH1F("h_ZpM_sw", "h_ZpM_sw", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    TH1F* h_zpM_ws    = new TH1F("h_ZpM_ws", "h_ZpM_ws", 50,Zpmass.Atof()-250,Zpmass.Atof()+250);
    
    TH1F* h_hPt       = new TH1F("h_higgsPt", "h_higgsPt", 60,0,1200);
    TH1F* h_a0Pt      = new TH1F("h_a0Pt", "h_A0Pt", 60,0,1200);
    TH1F* h_hPt_sw    = new TH1F("h_higgsPt_sw", "h_higgsPt_sw", 60,0,1200);
    TH1F* h_a0Pt_sw   = new TH1F("h_a0Pt_sw", "h_A0Pt_sw", 60,0,1200);
    TH1F* h_hPt_ws    = new TH1F("h_higgsPt_ws", "h_higgsPt_ws", 60,0,1200);
    TH1F* h_a0Pt_ws   = new TH1F("h_a0Pt_ws", "h_A0Pt_ws", 60,0,1200);
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
    TH1F *h_CISVV2_H[4], *h_CISVV2_A0[4];
    TH1F *h_TMVAweight_sw[2];
    TH1F *h_TMVAweight_ws[2];
    /*
    h_hM->SetLineColor(kBlue);
    h_zpM->SetLineColor(kBlue);
    h_hPt->SetLineColor(kBlue);
    h_hDeltaR->SetLineColor(kBlue);
    h_hDeltaPhi->SetLineColor(kBlue);
    h_hDeltaEta->SetLineColor(kBlue);
    h_hPtAs->SetLineColor(kBlue);
    h_hPtSD->SetLineColor(kBlue);
    h_zpDeltaR->SetLineColor(kBlue);
    h_zpDeltaPhi->SetLineColor(kBlue);
    h_zpDeltaR->SetLineColor(kBlue);
    h_zpDeltaEta->SetLineColor(kBlue);
    h_zpPtSD->SetLineColor(kBlue);
    h_zpPtAs->SetLineColor(kBlue);
    */
    for (int i=0;i<2;i++) {
        h_CISVV2_H[i]     = new TH1F(Form("h_CISVV2_Hsw_%d",i),Form("h_CISVV2_Hsw_%d",i),50,0.5,1);
        h_CISVV2_H[2+i]     = new TH1F(Form("h_CISVV2_Hws_%d",i),Form("h_CISVV2_Hws_%d",i),50,0.5,1);
        h_CISVV2_A0[i]    = new TH1F(Form("h_CISVV2_A0sw_%d",i),Form("h_CISVV2_A0sw_%d",i),50,0.5,1);
        h_CISVV2_A0[2+i]    = new TH1F(Form("h_CISVV2_A0ws_%d",i),Form("h_CISVV2_A0ws_%d",i),50,0.5,1);
        h_TMVAweight_ws[i]  = new TH1F(Form("h_TMVAweightWS_%d",i),Form("h_TMVAweightWS_%d",i),40,0.0,2);
        h_TMVAweight_sw[i]  = new TH1F(Form("h_TMVAweightSW_%d",i),Form("h_TMVAweightSW_%d",i),40,0.0,2);
        
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
    
    Int_t nPass[20]={0};
    int maxHIndexNum = 0, maxZpIndexNum = 0;
    
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        //if (jEntry>500) break; 
        nPass[0]++;
        float HT = data.GetFloat("HT");
        //h_HT->Fill(HT);
        h_allEvent->Fill(1);
        
        //0. has a good vertex
        int nGenJet =  data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[1]++;
        
        vector<bool>& vPassID_L = *(vector<bool>*) data.GetPtr("THINjetPassIDLoose");
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float *CISVV2 = data.GetPtrFloat("THINjetCISVV2");
        
        vector<vector<float>> HindexList, A0indexList, ZpindexList;
        vector<TLorentzVector*> genHA0Par;
        
        if (!isBG||doGenMatch) {
            TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
            for (int i=0;i<4;i++) genHA0Par.push_back((TLorentzVector*)genParP4->At(genPar[jEntry][i]));
        }
        
        // find good jet
        static vector<int> goodJet;
        goodJet.clear();
        static TLorentzVector* thisJet, *thatJet;
        for (int i=0;i<nGenJet;i++) {
            if (!vPassID_L[i]) continue;
            if (CISVV2[i]<CISVV2Cut) continue;
            thisJet = (TLorentzVector*)genjetP4->At(i);
            if (thisJet->Pt()<30) continue;
            if (thisJet->Eta()>2.4) continue;
            goodJet.push_back(i);
        }
        h_NgoodJets->Fill(goodJet.size());
        if (goodJet.size()<4) continue;
        nPass[2]++;
        TMVAinputInfo info(data);
        
        static vector<HA0JetInfo> fourJetList;
        static HA0JetInfo minX2, minWeightX2; 
        fourJetList.clear();
        static unsigned int hi, hj, ai, aj;
        /*
        for (hi=0;hi<goodJet.size();hi++) {
            for (hj=0;hj<hi;hj++) {
                for (ai=0;ai<goodJet.size();ai++) {
                    if ( ai==hi || ai==hj ) continue;
                    for (aj=0;aj<ai;aj++) {
                        if ( aj==hi || aj==hj ) continue;
                        fourJetList.push_back(HA0JetInfo(info,hi,hj,ai,aj,1));
                    }
                }
            }
        }
        */
        for (hi=0;hi<goodJet.size();hi++) {
            for (hj=0;hj<hi;hj++) {
                for (ai=0;ai<goodJet.size();ai++) {
                    if ( ai==hi || ai==hj ) continue;
                    for (aj=0;aj<ai;aj++) {
                        if ( aj==hi || aj==hj ) continue;
                        fourJetList.push_back(HA0JetInfo(info,goodJet[hi],goodJet[hj],goodJet[ai],goodJet[aj],1));
                    }
                }
            }
        }
        h_Ncom->Fill(fourJetList.size());
        
        if (fourJetList.size()==0) continue;
        nPass[3]++;

        sort(fourJetList.begin(),fourJetList.end(),sortListbyxSquare);
        minX2       = fourJetList[0];
        sort(fourJetList.begin(),fourJetList.end(),sortListbyWeightxSquare);
        minWeightX2 = fourJetList[0];
        
        // fill TH1F
        // Mass 
        h_hM->Fill(minX2.Mh());
        h_hM_sw->Fill(minX2.weightJets()->Mh());
        h_hM_ws->Fill(minWeightX2.weightJets()->Mh());
        h_a0M->Fill(minX2.MA0());
        h_a0M_sw->Fill(minX2.weightJets()->MA0());
        h_a0M_ws->Fill(minWeightX2.weightJets()->MA0());
        h_zpM->Fill(minX2.MZp());
        h_zpM_sw->Fill(minX2.weightJets()->MZp());
        h_zpM_ws->Fill(minWeightX2.weightJets()->MZp());
        // Pt
        h_hPt->Fill(minX2.Pth()); 
        h_hPt_sw->Fill(minX2.weightJets()->Pth()); 
        h_hPt_ws->Fill(minWeightX2.weightJets()->Pth());
        h_a0Pt->Fill(minX2.Pta0());
        h_a0Pt_sw->Fill(minX2.weightJets()->Pta0());
        h_a0Pt_ws->Fill(minWeightX2.weightJets()->Pta0());
        // Weight
        h_TMVAweight_sw[0]->Fill(minX2.getWeight(0)); 
        h_TMVAweight_sw[1]->Fill(minX2.getWeight(1)); 
        h_TMVAweight_ws[0]->Fill(minWeightX2.getWeight(0)); 
        h_TMVAweight_ws[1]->Fill(minWeightX2.getWeight(1));
        // CISVV2
        h_CISVV2_H[0]->Fill(CISVV2[minX2.getIndex(0)]);
        h_CISVV2_H[1]->Fill(CISVV2[minX2.getIndex(1)]);
        h_CISVV2_A0[0]->Fill(CISVV2[minX2.getIndex(2)]);
        h_CISVV2_A0[1]->Fill(CISVV2[minX2.getIndex(3)]);
        h_CISVV2_H[2]->Fill(CISVV2[minWeightX2.getIndex(0)]);
        h_CISVV2_H[3]->Fill(CISVV2[minWeightX2.getIndex(1)]);
        h_CISVV2_A0[2]->Fill(CISVV2[minWeightX2.getIndex(2)]);
        h_CISVV2_A0[3]->Fill(CISVV2[minWeightX2.getIndex(3)]);
        // Assymetry
        static TLorentzVector *Hjet, *Hjet_sw, *Hjet_ws, *A0jet, *A0jet_sw, *A0jet_ws;
        Hjet        = minX2.getHA0jet(0);
        Hjet_sw     = minX2.weightJets()->getHA0jet(0);
        Hjet_ws     = minWeightX2.weightJets()->getHA0jet(0);
        A0jet       = minX2.getHA0jet(1);
        A0jet_sw    = minX2.weightJets()->getHA0jet(1);
        A0jet_ws    = minWeightX2.weightJets()->getHA0jet(1);
        
        h_hPtAs->Fill(ptAssymetry(minX2.getHjet(0),minX2.getHjet(1)));
        h_hPtAs_sw->Fill(ptAssymetry(minX2.weightJets()->getHjet(0),minX2.weightJets()->getHjet(1)));
        h_hPtAs_ws->Fill(ptAssymetry(minWeightX2.weightJets()->getHjet(0),minWeightX2.weightJets()->getHjet(1)));
        h_a0PtAs->Fill(ptAssymetry(minX2.getA0jet(0),minX2.getA0jet(1)));
        h_a0PtAs_sw->Fill(ptAssymetry(minX2.weightJets()->getA0jet(0),minX2.weightJets()->getA0jet(1)));
        h_a0PtAs_ws->Fill(ptAssymetry(minWeightX2.weightJets()->getA0jet(0),minWeightX2.weightJets()->getA0jet(1)));
        h_zpPtAs->Fill(ptAssymetry(Hjet,A0jet));
        h_zpPtAs_sw->Fill(ptAssymetry(Hjet_sw,A0jet_sw));
        h_zpPtAs_ws->Fill(ptAssymetry(Hjet_ws,A0jet_sw));
        // SD 
        h_hPtSD->Fill(softDropAs(minX2.getHjet(0),minX2.getHjet(1)));
        h_hPtSD_sw->Fill(softDropAs(minX2.weightJets()->getHjet(0),minX2.weightJets()->getHjet(1)));
        h_hPtSD_ws->Fill(softDropAs(minWeightX2.weightJets()->getHjet(0),minWeightX2.weightJets()->getHjet(1)));
        h_a0PtSD->Fill(softDropAs(minX2.getA0jet(0),minX2.getA0jet(1)));
        h_a0PtSD_sw->Fill(softDropAs(minX2.weightJets()->getA0jet(0),minX2.weightJets()->getA0jet(1)));
        h_a0PtSD_ws->Fill(softDropAs(minWeightX2.weightJets()->getA0jet(0),minWeightX2.weightJets()->getA0jet(1)));
        h_zpPtSD->Fill(softDropAs(Hjet, A0jet));
        h_zpPtSD_sw->Fill(softDropAs(Hjet_sw,A0jet_sw));
        h_zpPtSD_ws->Fill(softDropAs(Hjet_ws,A0jet_ws));
        
        // Delta R, theta, phi
        // weight does not affect eta and phi
        h_hDeltaR_sw->Fill(minX2.getHjet(0)->DeltaR(*minX2.getHjet(1)));
        h_hDeltaR_ws->Fill(minWeightX2.getHjet(0)->DeltaR(*minWeightX2.getHjet(1)));
        h_a0DeltaR_sw->Fill(minX2.getA0jet(0)->DeltaR(*minX2.getHjet(1)));
        h_a0DeltaR_ws->Fill(minWeightX2.getA0jet(0)->DeltaR(*minWeightX2.getHjet(1)));
        h_zpDeltaR_sw->Fill(Hjet_sw->DeltaR(*A0jet_sw));
        h_zpDeltaR_ws->Fill(Hjet_ws->DeltaR(*A0jet_ws));

        h_hDeltaEta_sw->Fill(abs(minX2.getHjet(0)->Eta() - minX2.getHjet(1)->Eta())); 
        h_hDeltaEta_ws->Fill(abs(minWeightX2.getHjet(0)->Eta() - minWeightX2.getHjet(1)->Eta())); 
        h_a0DeltaEta_sw->Fill(abs(minX2.getA0jet(0)->Eta() - minX2.getA0jet(1)->Eta())); 
        h_a0DeltaEta_ws->Fill(abs(minWeightX2.getA0jet(0)->Eta() - minWeightX2.getA0jet(1)->Eta())); 
        h_zpDeltaEta_sw->Fill(abs(Hjet_sw->Eta()-A0jet_sw->Eta()));
        h_zpDeltaEta_ws->Fill(abs(Hjet_ws->Eta()-A0jet_ws->Eta()));

        h_hDeltaPhi_sw->Fill(caldePhi(minX2.getHjet(0)->Phi(),minX2.getHjet(1)->Phi()));
        h_hDeltaPhi_ws->Fill(caldePhi(minWeightX2.getHjet(0)->Phi(),minWeightX2.getHjet(1)->Phi()));
        h_a0DeltaPhi_sw->Fill(caldePhi(minX2.getA0jet(0)->Phi(),minX2.getA0jet(1)->Phi()));
        h_a0DeltaPhi_ws->Fill(caldePhi(minWeightX2.getA0jet(0)->Phi(),minWeightX2.getA0jet(1)->Phi()));
        h_zpDeltaPhi_sw->Fill(caldePhi(Hjet_sw->Phi(),A0jet_sw->Phi()));
        h_zpDeltaPhi_ws->Fill(caldePhi(Hjet_ws->Phi(),A0jet_ws->Phi()));
        // 
        //h_TMVAweight_ws[0]->Fill(minX2.weightJets()->getWeight(0)); 
        //h_TMVAweight_ws[1]->Fill(minX2.weightJets()->getWeight(1));
        //cout << minX2.getWeight(0) << " " << minX2.getWeight(1)  << "| \t" << minX2.weightJets()->getWeight(0) << " " << minX2.weightJets()->getWeight(1) << endl;
        //cout << minX2.getWeight(0) << " " << minX2.getWeight(1)  << "| \t" << minWeightX2.getWeight(0) << " " << minWeightX2.getWeight(1) << endl;
        
        
        
       delete Hjet;
       delete Hjet_sw;
       delete Hjet_ws;
       delete A0jet;
       delete A0jet_sw;
       delete A0jet_ws;
        
    } // end of loop over entries
    /*
    h_zpPtSD->SetMaximum(1000);
    h_zpPtAs_ori->SetMaximum(150); 
    */
    float nTotal = data.GetEntriesFast();
    std::cout << "nTotal    = " << nTotal << std::endl;
    
    for(int i=0;i<20;i++) if(nPass[i]>0) std::cout << "nPass[" << i << "] = " << nPass[i] << std::endl;
    efferr(nPass[3],nTotal);
    TLegend* leg = new TLegend(0.68,0.8,0.9,0.9);
    //leg->AddEntry(h_hM,"TMVA Reweight");
    //leg->AddEntry(h_hM_ori,"origin");
    if (!isBG) 
    {
        string pdfName;
        if (doGenMatch) pdfName = Form("anGenMatchJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        else if (doBTagCut) pdfName = Form("anRecoBTagJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        else pdfName = Form("anRecoJet_bb2HDM_MZp%s_MA0%s.pdf",Zpmass.Data(),A0mass.Data());
        c1->Print((pdfName+"[").data());
        h_NgoodJets->Draw("hist");
        c1->Print(pdfName.data());
        h_Ncom->Draw("hist");
        c1->Print(pdfName.data());
        draw(h_hM,h_hM_sw,h_hM_ws);
        c1->Print(pdfName.data());
        h_a0M->Draw("hist");
        c1->Print(pdfName.data());
        h_a0M_sw->Draw("hist");
        c1->Print(pdfName.data());
        h_a0M_ws->Draw("hist");
        c1->Print(pdfName.data());
        draw(h_a0M,h_a0M_sw,h_a0M_ws);
        c1->Print(pdfName.data());
        draw(h_zpM,h_zpM_sw,h_zpM_ws);
        c1->Print(pdfName.data());
        h_a0Pt->SetTitle("h_HA0rawPt_comparison");
        draw(h_a0Pt,h_hPt);
        c1->Print(pdfName.data());
        h_a0Pt_sw->SetTitle("h_HA0PtSW_comparison");
        draw(h_a0Pt_sw,h_hPt_sw);
        c1->Print(pdfName.data());
        h_a0Pt_ws->SetTitle("h_HA0PtWS_comparison");
        draw(h_a0Pt_ws,h_hPt_ws);
        c1->Print(pdfName.data());
        h_hPt->SetTitle("h_hPt_raw_WS_SW");
        draw(h_hPt,h_hPt_ws,h_hPt_sw);
        c1->Print(pdfName.data());
        h_hPt->SetTitle("h_a0Pt_raw_WS_SW");
        draw(h_a0Pt,h_a0Pt_ws,h_a0Pt_sw);
        c1->Print(pdfName.data());
        h_CISVV2_H[0]->SetTitle("h_hCISVV2_SW");
        draw(h_CISVV2_H[0],h_CISVV2_H[1]);
        c1->Print(pdfName.data());
        h_CISVV2_H[2]->SetTitle("h_hCISVV2_WS");
        draw(h_CISVV2_H[2],h_CISVV2_H[3]);
        c1->Print(pdfName.data());
        h_CISVV2_A0[0]->SetTitle("h_a0CISVV2_SW");
        draw(h_CISVV2_A0[0],h_CISVV2_A0[1]);
        c1->Print(pdfName.data());
        h_CISVV2_A0[2]->SetTitle("h_a0CISVV2_WS");
        draw(h_CISVV2_A0[2],h_CISVV2_A0[3]);
        c1->Print(pdfName.data());
        h_TMVAweight_sw[0]->SetTitle("h_hWeight_SW");
        draw(h_TMVAweight_sw[0],h_TMVAweight_sw[1]);
        c1->Print(pdfName.data());
        h_TMVAweight_sw[0]->SetTitle("h_hWeight_WS");
        draw(h_TMVAweight_ws[0],h_TMVAweight_ws[1]);
        c1->Print(pdfName.data());
        draw(h_hPtAs,h_hPtAs_sw,h_hPtAs_ws);
        c1->Print(pdfName.data());
        draw(h_a0PtAs,h_a0PtAs_sw,h_a0PtAs_ws);
        c1->Print(pdfName.data());
        draw(h_hPtSD,h_hPtSD_sw,h_hPtSD_ws);
        c1->Print(pdfName.data());
        draw(h_a0PtSD,h_a0PtSD_sw,h_a0PtSD_ws);
        c1->Print(pdfName.data());
        draw(h_hDeltaR_sw,h_hDeltaR_ws);
        c1->Print(pdfName.data());
        draw(h_a0DeltaR_sw,h_a0DeltaR_ws);
        c1->Print(pdfName.data());
        draw(h_zpDeltaR_sw,h_zpDeltaR_ws);
        c1->Print(pdfName.data());
        draw(h_hDeltaEta_sw,h_hDeltaEta_ws);
        c1->Print(pdfName.data());
        draw(h_a0DeltaEta_sw,h_a0DeltaEta_ws);
        c1->Print(pdfName.data());
        draw(h_zpDeltaEta_sw,h_zpDeltaEta_ws);
        c1->Print(pdfName.data());
        draw(h_hDeltaPhi_sw,h_hDeltaPhi_ws);
        c1->Print(pdfName.data());
        draw(h_a0DeltaPhi_sw,h_a0DeltaPhi_ws);
        c1->Print(pdfName.data());
        draw(h_zpDeltaPhi_sw,h_zpDeltaPhi_ws);
        c1->Print(pdfName.data());

        c1->Print((pdfName+"]").data());
    }
    /*
    string fileName; 
    if (isBG) fileName = Form("QCDbg2HDMbb_%d.root",w);
    else if (doBTagCut) fileName = Form("sigRootFile/bb2HDM_recoBTag_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    else fileName = Form("sigRootFile/bb2HDM_reco_MZp%s_MA0%s.root",Zpmass.Data(),A0mass.Data());
    if (!saveEffi&&false) {
        TFile* outputFile = new TFile(fileName.data(),"recreate");
        h_HT->Write();
        h_zpM->Write();
        h_hM->Write();
        h_a0M->Write();
        h_hPt->Write();
        h_a0Pt->Write();
        h_hPtAs->Write();
        h_a0PtAs->Write();
        h_zpPtAs->Write();
        h_hPtSD->Write();
        h_a0PtSD->Write();
        h_zpPtSD->Write();
        h_hDeltaR->Write();
        h_a0DeltaR->Write();
        h_zpDeltaR->Write();
        h_hDeltaEta->Write();
        h_a0DeltaEta->Write();
        h_zpDeltaEta->Write();
        h_hDeltaPhi->Write();
        h_a0DeltaPhi->Write();
        h_zpDeltaPhi->Write();
        h_CISVV2_H[0]->Write();
        h_CISVV2_H[1]->Write();
        h_CISVV2_A0[0]->Write();
        h_CISVV2_A0[1]->Write();
        h_hNcandi->Write();
        h_a0Ncandi->Write();
        h_zpNcandi->Write();
        
        outputFile->Close();
    }
    string fileNPassName;
    if (doBTagCut) fileNPassName = Form("../njetsnPass/effi_Zpmass%s_A0mass%s_recoBTag4jets.txt",Zpmass.Data(),A0mass.Data());
    else fileNPassName = Form("../njetsnPass/effi_Zpmass%s_A0mass%s_reco4jets.txt",Zpmass.Data(),A0mass.Data());
    savenPass(nPass,fileNPassName);
    if (saveEffi) {
        string effifilename;
        if (doBTagCut) effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_recoBTag4jets.txt",Zpmass.Data(),A0mass.Data());
        else effifilename = Form("../njetsEffi/effi_Zpmass%s_A0mass%s_reco4jets.txt",Zpmass.Data(),A0mass.Data());
        //string nPassfilename = Form("../njetsnPass/effi_Zpmass%s_A0mass%s_4jets.txt",Zpmass.Data(),A0mass.Data());
        fstream fp;
        fp.open(effifilename.data(), ios::out);
        fp << (float)nPass[4]/nTotal << endl;
        fp.close();
        
    }*/
    delete leg;
}
