#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include "string"
#include "setNCUStyle.C"
#include "pfRatio.h"

#define basePtEtaCut true
using namespace std;

bool sortfourJets(fourJetInfo a, fourJetInfo b) {return a.ChiSquare<b.ChiSquare;}
void pfRatio(string inputFile="../../QCDtestBGrootfile/NCUGlobalTuples_76.root") {
    
    setNCUStyle(true);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0001111101.);
    int w = 0;
    
    const float CISVV2CUT_L = 0.5426;
    const float CISVV2CUT_M = 0.8484;
    const float CISVV2CUT_T = 0.9535;
    
    TFile *f = TFile::Open(Form("pfRatio_%d.root",w),"RECREATE");

    TCanvas* c1 = new TCanvas("c1","",600,600);
    const int nCom = 20; 
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcombinationJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcombinationJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcombinationJets", nCom,-0.5,nCom-0.5);
    TH1F* h_A0M_noCISVV2 = new TH1F("h_A0M_noCISVV2","h_A0M_noCISVV2",100,-50,50);
    TH1F* h_A0M_CISVV2 = new TH1F("h_A0M_CISVV2","h_A0M_CISVV2",100,-50,50);
    TH1F* h_pfA0M[2];
    h_pfA0M[0] = new TH1F("h_pfA0M_f","h_pfA0M_f",150,225,275);
    h_pfA0M[1] = new TH1F("h_pfA0M_p","h_pfA0M_p",150,225,275);
    TH1F* h_pfCISVV2[2];
    h_pfCISVV2[0] = new TH1F("h_pfCISVV2_f","h_pfCISVV2_f",150,-75,75);
    h_pfCISVV2[1] = new TH1F("h_pfCISVV2_p","h_pfCISVV2_p",150,-75,75);
    TH1F* h_pfCISVV2_ratio = new TH1F("h_pfCISVV2_ratio","h_pfCISVV2_ratio",150,-75,75);

    vector<vector<float>> HindexList, A0indexList;
    vector<fourJetInfo> ZpindexList;
    fourJetInfo fourJets;
    Int_t nPass[20]={0};
    
    TreeReader data(inputFile.data());
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++) {
        if (jEntry %2000 == 0) fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        data.GetEntry(jEntry);
        nPass[0]++;
        int nGenJet =  data.GetInt("THINnJet");
        if(nGenJet<4) continue;
        nPass[1]++;
        //static vector<bool> isHpassCISVV2;
        static bool isHpassCISVV2 = false;
        static bool isA0passCISVV2 = false;
        TClonesArray* genjetP4 =  (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float *CISVV2 = data.GetPtrFloat("THINjetCISVV2");
        HindexList.clear();
        ZpindexList.clear();

        for(int ij=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if (thisJet->Pt()<30)continue;
            if (fabs(thisJet->Eta())>2.4)continue;
            //if (CISVV2[ij]<CISVV2CUT_L) continue;
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if(thatJet->Pt()<30)continue;
                if(fabs(thatJet->Eta())>2.4)continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (diJetM>140||diJetM<110) continue;
                if (CISVV2[jj]<CISVV2CUT_L||CISVV2[ij]<CISVV2CUT_L) isHpassCISVV2 = false;
                else isHpassCISVV2 = true;
                HindexList.push_back({(float)ij,(float)jj,(float)(*thisJet+*thatJet).Pt(),(float)isHpassCISVV2});
            }
        } // end of outer loop jet
        sort(HindexList.begin(),HindexList.end(),sortListbyPt);
        h_hNcandi->Fill(HindexList.size());
        if (HindexList.size()==0) continue;
        nPass[2]++;
    
        // choose A0 and apply Zp mass cut
        for(int i=0;i<HindexList.size();i++) {
            for(int ij=0; ij < nGenJet; ij++) {
                // skip repeat jet
                if (ij==(int)HindexList[i][0]) continue;
                if (ij==(int)HindexList[i][1]) continue;
                
                TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
                if(basePtEtaCut && thisJet->Pt()<30)continue;
                if(basePtEtaCut && fabs(thisJet->Eta())>2.4)continue;
                for (int jj=0;jj<ij;jj++) {
                    // skip repeat jet
                    if (jj==(int)HindexList[i][1]) continue;
                    
                    TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                    if(basePtEtaCut && thatJet->Pt()<30) continue;
                    if(basePtEtaCut && fabs(thatJet->Eta())>2.4) continue;
                    TLorentzVector* HJet1 = (TLorentzVector*)genjetP4->At((int)HindexList[i][0]);
                    TLorentzVector* HJet2 = (TLorentzVector*)genjetP4->At((int)HindexList[i][1]);
                    //fourJets = fourJetInfo(HJet1,HJet2,thisJet,thatJet);
                    fourJets = fourJetInfo(data,(int)HindexList[i][0],(int)HindexList[i][1],ij,jj);
                    if (fourJets.MZp<900 || fourJets.MZp>1100) continue;
                    ZpindexList.push_back(fourJets); 
                }
            } // end of A0 loop jet
        }// end of H index loop
        h_zpNcandi->Fill(ZpindexList.size());
        if (ZpindexList.size()==0) continue;
        sort(ZpindexList.begin(),ZpindexList.end(),sortfourJets);
        nPass[3]++;
        isA0passCISVV2 = ZpindexList[0].minCISVV2[1]>CISVV2CUT_L;
        h_pfCISVV2[isA0passCISVV2]->Fill(ZpindexList[0].MA0-300);
        h_pfA0M[isA0passCISVV2]->Fill(ZpindexList[0].MA0);
        //if(jEntry>=2000) break;
    } // end of event loop
    
    h_pfCISVV2_ratio->Divide(h_pfCISVV2[1],h_pfCISVV2[0]);
    string pdfName = "tes.pdf";
    
    c1->Print((pdfName+"[").data());
    h_pfCISVV2_ratio->Draw("hist");
    c1->Print(pdfName.data());
    h_pfCISVV2[0]->Draw("hist");
    c1->Print(pdfName.data());
    h_pfCISVV2[1]->Draw("hist");
    c1->Print(pdfName.data());
    h_pfA0M[0]->Draw("hist");
    c1->Print(pdfName.data());
    h_pfA0M[1]->Draw("hist");
    c1->Print(pdfName.data());
    h_hNcandi->Draw("hist text 0");
    c1->Print(pdfName.data());
    h_zpNcandi->Draw("hist text 0");
    c1->Print(pdfName.data());

    c1->Print((pdfName+"]").data());
}
