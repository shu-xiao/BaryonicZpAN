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
#define savePDF true

#define isSig 0
#if isSig
    #define INPUTFILE "../2HDMfullSimFile/2HDM_MZp1000_MA0300.root"
#else
    #define INPUTFILE "../../QCDtestBGrootfile/NCUGlobalTuples_76.root"
#endif
using namespace std;

bool sortfourJets(fourJetInfo a, fourJetInfo b) {return a.ChiSquare<b.ChiSquare;}
void pfRatio(int w=0,string inputFile=INPUTFILE) {
    
    setNCUStyle(true);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0001112101.);
    gStyle->SetOptTitle(1);
    
    const float CISVV2CUT_L = 0.5426;
    const float CISVV2CUT_M = 0.8484;
    const float CISVV2CUT_T = 0.9535;
    
    bool isBG = false;
    
    TString Zpmass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*_MZp}; test=${test1%%_MA0*.root}; echo \"${test}\"",inputFile.data()));
    TString A0mass=gSystem->GetFromPipe(Form("file=%s; test1=${file##*MA0}; test=${test1%%.root}; echo \"${test}\"",inputFile.data()));
    // set default value
    if (Zpmass.Atof()>2500||Zpmass.Atof()<500 || A0mass.Atof()>800 || A0mass.Atof()<300) {
        Zpmass = "1000";
        A0mass = "300";
        isBG = true;
    }
    TFile *f = TFile::Open(Form("pfRatio_%d.root",w),"RECREATE");

    TCanvas* c1 = new TCanvas("c1","",600,600);
    const int nCom = 20; 
    TH1F* h_hNcandi = new TH1F("h_hNcandi", "h_higgs_NcombinationJets", 10,-0.5,9.5);
    TH1F* h_a0Ncandi = new TH1F("h_a0Ncandi", "h_A0_NcombinationJets", 10,-0.5,9.5);
    TH1F* h_zpNcandi = new TH1F("h_zpNcandi", "h_Zp_NcombinationJets", nCom,-0.5,nCom-0.5);
    TH1F* h_A0M_noCISVV2 = new TH1F("h_A0M_noCISVV2","h_A0M_noCISVV2",30,-75,75);
    TH1F* h_A0M_CISVV2 = new TH1F("h_A0M_CISVV2","h_A0M_CISVV2",30,-75,75);
    TH1F* h_pfA0M[2];
    h_pfA0M[0] = new TH1F("h_pfA0M_F","h_pfA0M_F",30,A0mass.Atof()-75,A0mass.Atof()+75);
    h_pfA0M[1] = new TH1F("h_pfA0M_P","h_pfA0M_P",30,A0mass.Atof()-75,A0mass.Atof()+75);
    TH1F* h_sigA0M[9];
    for (int i=0;i<9;i++) h_sigA0M[i] = new TH1F(Form("h_sigA0M_%d",i),Form("h_sigA0M_%d",i),30,A0mass.Atof()-75,A0mass.Atof()+75);
    TH1F* h_pfCISVV2[2];
    h_pfCISVV2[0] = new TH1F("h_pfCISVV2_F","h_pfCISVV2_F",30,-75,75);
    h_pfCISVV2[1] = new TH1F("h_pfCISVV2_P","h_pfCISVV2_P",30,-75,75);
    TH1F* h_pfCISVV2_ratio = new TH1F("h_pfCISVV2_ratio","h_pfCISVV2_ratio",30,-75,75);

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

        // find higgs
        for(int ij=0; ij < nGenJet; ij++) {
            TLorentzVector* thisJet = (TLorentzVector*)genjetP4->At(ij);
            if (thisJet->Pt()<30)continue;
            if (fabs(thisJet->Eta())>2.4)continue;
            if (CISVV2[ij]<CISVV2CUT_L) continue;
            for (int jj=0;jj<ij;jj++) {
                TLorentzVector* thatJet = (TLorentzVector*)genjetP4->At(jj);
                if(thatJet->Pt()<30)continue;
                if(fabs(thatJet->Eta())>2.4)continue;
                float diJetM = (*thisJet+*thatJet).M();
                if (CISVV2[jj]<CISVV2CUT_L) continue;
                if (diJetM>140||diJetM<110) continue;
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
                    if (fourJets.MZp<Zpmass.Atof()-100 || fourJets.MZp>Zpmass.Atof()+100) continue;
                    ZpindexList.push_back(fourJets); 
                }
            } // end of A0 loop jet
        }// end of H index loop
        h_zpNcandi->Fill(ZpindexList.size());
        if (ZpindexList.size()==0) continue;
        sort(ZpindexList.begin(),ZpindexList.end(),sortfourJets);
        nPass[3]++;
        
        // take pass or fail event
        isA0passCISVV2 = ZpindexList[0].minCISVV2[1]>CISVV2CUT_L; // choose highest possible A0 jet min CISVV2
        h_pfCISVV2[isA0passCISVV2]->Fill(ZpindexList[0].MA0 - A0mass.Atof());
        h_pfA0M[isA0passCISVV2]->Fill(ZpindexList[0].MA0);

        // observe the higgs mass CISVV2 cut 
        static int i;
        // LL MM TT
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].CISVV2_HA0[0]>CISVV2CUT_T&&ZpindexList[i].CISVV2_HA0[1]>CISVV2CUT_T){
                h_sigA0M[0]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].CISVV2_HA0[0]>CISVV2CUT_M&&ZpindexList[i].CISVV2_HA0[1]>CISVV2CUT_M){
                h_sigA0M[1]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].CISVV2_HA0[0]>CISVV2CUT_L&&ZpindexList[i].CISVV2_HA0[1]>CISVV2CUT_L){
                h_sigA0M[2]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        // min CISVV2
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].minCISVV2[0]>CISVV2CUT_T){
                h_sigA0M[3]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].minCISVV2[0]>CISVV2CUT_M){
                h_sigA0M[4]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].minCISVV2[0]>CISVV2CUT_L){
                h_sigA0M[5]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        // maxPtCISVV2
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].maxPtCISVV2[0]>CISVV2CUT_T){
                h_sigA0M[6]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].maxPtCISVV2[0]>CISVV2CUT_M){
                h_sigA0M[7]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        for (i=0;i<ZpindexList.size();i++) {
            if (ZpindexList[i].maxPtCISVV2[0]>CISVV2CUT_L){
                h_sigA0M[8]->Fill(ZpindexList[i].MA0);
                break;
            }
        }
        //if(jEntry>=2000) break;
    } // end of event loop
    
    h_pfCISVV2_ratio->Divide(h_pfCISVV2[1],h_pfCISVV2[0]);
    string pdfName;
    if (!isBG) pdfName = "tes_sig.pdf";
    else pdfName = "tes.pdf";
    int colorInd[3] = {kBlue+3,kBlue,kBlue-3}; 
    string title[3] = {"TT MM LL","Min CISVV2","Max Pt CISVV2"};
    if (savePDF) {
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
        for (int i=0;i<3;i++) {
            c1->Clear();
            for (int j=0;j<3;j++) h_sigA0M[3*i+j]->SetLineColor(colorInd[2-j]);
            auto leg = TLegend(0.7,0.5,0.95,0.7);
            leg.AddEntry(h_sigA0M[3*i+2],(i)?Form("T, RMS=%.1f",h_sigA0M[3*i+2]->GetRMS()):Form("TT, RMS=%.1f",h_sigA0M[3*i+2]->GetRMS()));
            leg.AddEntry(h_sigA0M[3*i+1],(i)?Form("M, RMS=%.1f",h_sigA0M[3*i+1]->GetRMS()):Form("MM, RMS=%.1f",h_sigA0M[3*i+1]->GetRMS()));
            leg.AddEntry(h_sigA0M[3*i],(i)?Form("L, RMS=%.1f",h_sigA0M[3*i]->GetRMS()):Form("LL, RMS=%.1f",h_sigA0M[3*i]->GetRMS()));
            h_sigA0M[3*i]->SetTitle(title[i].data());
            h_sigA0M[3*i+1]->SetTitle(title[i].data());
            h_sigA0M[3*i+2]->SetTitle(title[i].data());
            h_sigA0M[3*i+2]->Draw("hist");
            h_sigA0M[3*i+1]->Draw("histsame");
            h_sigA0M[3*i]->Draw("histsame");
            leg.Draw();
            c1->Print(pdfName.data());
            
        }
        h_hNcandi->Draw("hist text 0");
        c1->Print(pdfName.data());
        h_zpNcandi->Draw("hist text 0");
        c1->Print(pdfName.data());

        c1->Print((pdfName+"]").data());
    }
    f->Write(); 
    f->Close();
}
